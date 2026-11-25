import test from 'node:test';
import assert from 'node:assert/strict';
import { layoutCaption } from '../render/overlay.js';

const measure = (s) => s.length * 10;

test('wraps long captions', () => {
  const text = 'the quick brown fox jumps over';
  const lines = layoutCaption(text, 300, measure);
  assert.ok(lines.length > 1);
  for (const l of lines) assert.ok(measure(l) <= 240, `"${l}" is ${measure(l)}px`);
  assert.equal(lines.join(' '), text);
});

test('single short line unchanged', () => {
  assert.deepEqual(layoutCaption('Hi', 1000, measure), ['Hi']);
});

test('a single word longer than the width stays on its own line', () => {
  assert.deepEqual(layoutCaption('a supercalifragilistic b', 100, measure), ['a', 'supercalifragilistic', 'b']);
});

test('explicit line breaks are kept', () => {
  assert.deepEqual(layoutCaption('one\ntwo', 1000, measure), ['one', 'two']);
});
