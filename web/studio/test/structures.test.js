import test from 'node:test';
import assert from 'node:assert/strict';
import { parse } from '../parse.js';
import { loadStructures, clearStructureCache } from '../structures.js';
import { chainPdb } from './fixtures.js';

const scene = (text) => {
  const r = parse(text);
  assert.deepEqual(r.errors, []);
  return r.scene;
};
const okJson = (body) => ({ ok: true, status: 200, json: async () => body });

test('loads via /api/analyze with chain', async () => {
  clearStructureCache();
  const urls = [];
  const fetchFn = async (url) => { urls.push(url); return okJson({ pdb_data: chainPdb('A', 12) }); };
  const r = await loadStructures(scene('actor sep = protein from=7NJ1 chain=A\nactor g = protein'), fetchFn);
  assert.deepEqual(urls, ['/api/analyze?pdb=7NJ1&chain=A']);
  assert.equal(r.beads.get('sep').beads.length, 3);
  assert.equal(r.beads.has('g'), false);
  assert.deepEqual(r.warnings, []);
});

test('falls back on 404 / missing chain', async () => {
  clearStructureCache();
  const r = await loadStructures(scene('\nactor p = protein from=9ZZZ chain=A'), async () => ({ ok: false, status: 404, json: async () => ({ error: 'not found' }) }));
  assert.equal(r.beads.get('p'), null);
  assert.equal(r.warnings.length, 1);
  assert.equal(r.warnings[0].line, 2);
  assert.match(r.warnings[0].message, /9ZZZ chain A could not be loaded .*showing a generic shape/);

  clearStructureCache();
  const big = await loadStructures(scene('actor p = protein from=6WG3 chain=A'), async () => okJson({ too_large: true }));
  assert.equal(big.beads.get('p'), null);
  assert.match(big.warnings[0].message, /showing a generic shape/);

  clearStructureCache();
  const none = await loadStructures(scene('actor p = protein from=7NJ1 chain=Q'), async () => okJson({ pdb_data: chainPdb('A', 4) }));
  assert.equal(none.beads.get('p'), null);
  assert.match(none.warnings[0].message, /no CA atoms/);

  clearStructureCache();
  const thrown = await loadStructures(scene('actor p = protein from=7NJ1'), async () => { throw new TypeError('Failed to fetch'); });
  assert.equal(thrown.beads.get('p'), null);
  assert.match(thrown.warnings[0].message, /Failed to fetch/);
});

test('multi-chain concatenates', async () => {
  clearStructureCache();
  const urls = [];
  const fetchFn = async (url) => {
    urls.push(url);
    const chain = new URL(url, 'http://x').searchParams.get('chain');
    return okJson({ pdb_data: chainPdb(chain, 8, chain === 'B' ? 40 : 0) });
  };
  const r = await loadStructures(scene('actor pol = protein from=1I6H chain=A,B'), fetchFn);
  assert.equal(urls.length, 2);
  assert.equal(r.beads.get('pol').beads.length, 4);
});

test('cache avoids refetching the same structure', async () => {
  clearStructureCache();
  let calls = 0;
  const fetchFn = async () => { calls++; return okJson({ pdb_data: chainPdb('A', 8) }); };
  await loadStructures(scene('actor a = protein from=7NJ1 chain=A\nactor b = protein from=7NJ1 chain=A'), fetchFn);
  await loadStructures(scene('actor a = protein from=7NJ1 chain=A'), fetchFn);
  assert.equal(calls, 1);
});
