import test from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync, readdirSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { parse } from '../parse.js';
import { compile } from '../timeline.js';

const DIR = fileURLToPath(new URL('../examples/', import.meta.url));
const EXPECTED = ['nondisjunction-meiosis-1.txt', 'transcription-basic.txt'];

test('the bundled examples are exactly the two the menu offers', () => {
  assert.deepEqual(readdirSync(DIR).filter((f) => f.endsWith('.txt')).sort(), EXPECTED);
  const html = readFileSync(fileURLToPath(new URL('../../studio.html', import.meta.url)), 'utf8');
  for (const f of EXPECTED) assert.ok(html.includes(`value="${f.replace('.txt', '')}"`), f);
});

for (const file of EXPECTED) {
  test(`${file} parses and compiles cleanly`, () => {
    const text = readFileSync(DIR + file, 'utf8');
    const p = parse(text);
    assert.deepEqual(p.errors, [], 'parse errors');
    const c = compile(p.scene);
    assert.deepEqual(c.errors, [], 'compile errors');
    assert.deepEqual([...p.warnings, ...c.warnings], [], 'warnings');
    const d = c.timeline.duration;
    assert.ok(d >= 30 && d <= 90, `duration ${d}s`);
    assert.equal(c.timeline.style, 'dark');
    assert.equal(c.timeline.note, 'Stylised; not to scale');
    assert.ok(c.timeline.creditsText, 'credits list real structures');
    assert.ok(c.timeline.captions.length >= 5, 'captions narrate each step');
    for (const a of p.scene.actors.filter((x) => x.type === 'protein' && x.params.from)) {
      assert.ok(String(a.params.chain).split(',').length <= 2, `${a.name}: at most two chains`);
    }
  });
}
