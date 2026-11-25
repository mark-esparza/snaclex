import test from 'node:test';
import assert from 'node:assert/strict';
import { caCoords, coarseBeads } from '../pdbcoarse.js';
import { pdbLine } from './fixtures.js';

const near = (a, b, eps = 1e-9) => assert.ok(Math.abs(a - b) < eps, `${a} ≉ ${b}`);

test('caCoords reads CA only, in nm, for one chain', () => {
  const pdb = [
    pdbLine('ATOM', 1, 'N', 'ALA', 'A', 1, 10, 0, 0),
    pdbLine('ATOM', 2, 'CA', 'ALA', 'A', 1, 11, 2, 3),
    pdbLine('ATOM', 3, 'C', 'ALA', 'A', 1, 12, 0, 0),
    pdbLine('ATOM', 4, 'CA', 'GLY', 'A', 2, 14, 0, 0),
    pdbLine('ATOM', 5, 'CA', 'SER', 'A', 3, 17, 0, 0),
    pdbLine('ATOM', 6, 'CA', 'SER', 'B', 1, 99, 0, 0),
    pdbLine('HETATM', 7, 'CA', ' CA', 'A', 900, 50, 50, 50),
  ].join('\n');
  const a = caCoords(pdb, 'A');
  assert.equal(a.length, 3);
  near(a[0][0], 1.1); near(a[0][1], 0.2); near(a[0][2], 0.3);
  assert.equal(caCoords(pdb, 'B').length, 1);
  assert.equal(caCoords(pdb, 'Z').length, 0);
});

test('caCoords reads only the first model and the first altloc', () => {
  const pdb = [
    'MODEL        1',
    pdbLine('ATOM', 1, 'CA', 'ALA', 'A', 1, 10, 0, 0, 'A'),
    pdbLine('ATOM', 2, 'CA', 'ALA', 'A', 1, 10.5, 0, 0, 'B'),
    'ENDMDL',
    'MODEL        2',
    pdbLine('ATOM', 1, 'CA', 'ALA', 'A', 1, 20, 0, 0),
    'ENDMDL',
  ].join('\n');
  assert.deepEqual(caCoords(pdb, 'A'), [[1, 0, 0]]);
});

test('coarseBeads centres and sizes', () => {
  const coords = Array.from({ length: 8 }, (_, i) => [i, 0, 0]);
  const b = coarseBeads(coords, 4);
  assert.equal(b.beads.length, 2);
  near(b.beads[0][0] + b.beads[1][0], 0);
  assert.equal(b.beads[0][3], 0.9);
  assert.ok(b.radius > 0.9);
  near(b.radius, 2 + 0.9);
});

test('coarseBeads keeps a short tail as its own bead', () => {
  const b = coarseBeads(Array.from({ length: 5 }, (_, i) => [i, 0, 0]), 4);
  assert.equal(b.beads.length, 2);
});
