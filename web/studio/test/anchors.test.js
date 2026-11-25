import test from 'node:test';
import assert from 'node:assert/strict';
import { rotate, actorExtent, anchorPoint, helixFrame } from '../anchors.js';

const close = (a, b, eps = 1e-9) => {
  assert.equal(a.length, b.length);
  a.forEach((v, i) => assert.ok(Math.abs(v - b[i]) < eps, `${a} ≉ ${b}`));
};
const state = (over = {}) => ({ pos: [0, 0, 0], rot: [0, 0, 0], scale: 1, p: {}, ...over });

test('rotate applies XYZ Euler order', () => {
  close(rotate([17, 0, 0], [0, 90, 0]), [0, 0, -17]);
  close(rotate([1, 0, 0], [0, 0, 90]), [0, 1, 0]);
  close(rotate([0, 1, 0], [90, 0, 0]), [0, 0, 1]);
  close(rotate([1, 2, 3], [0, 0, 0]), [1, 2, 3]);
});

test('helix anchor', () => {
  const decl = { type: 'helix', params: { turns: 10, curve: 0 } };
  close(anchorPoint(decl, state(), { ref: 'dna', u: 0 }), [-17, 0, 0]);
  close(anchorPoint(decl, state(), { ref: 'dna', u: 1 }), [17, 0, 0]);
  close(anchorPoint(decl, state({ rot: [0, 90, 0] }), { ref: 'dna', u: 1 }), [0, 0, -17]);
  close(anchorPoint(decl, state({ pos: [1, 2, 3], scale: 2 }), { ref: 'dna', u: 1 }), [35, 2, 3]);
  close(anchorPoint(decl, state(), { ref: 'dna', part: 'end' }), [17, 0, 0]);
  close(anchorPoint(decl, state(), { ref: 'dna' }), [0, 0, 0]);
});

test('helix anchor follows a cut', () => {
  const decl = { type: 'helix', params: { turns: 10, curve: 0 } };
  const s = state({ p: { cut: 1, cutAt: 0.5, gap: 4 } });
  close(anchorPoint(decl, s, { ref: 'dna', u: 1 }), [19, 0, 0]);
  close(anchorPoint(decl, s, { ref: 'dna', u: 0 }), [-19, 0, 0]);
});

test('curved helix keeps its midpoint at the origin', () => {
  const f = helixFrame(0.5, 10, 0.5);
  close(f.point, [0, 0, 0]);
  close(f.tangent, [1, 0, 0]);
  const end = helixFrame(1, 10, 0.5).point;
  assert.ok(end[1] > 0, 'ends bend upward');
});

test('chromosome parts', () => {
  const decl = { type: 'chromosome', params: { length: 12 } };
  close(anchorPoint(decl, state({ p: { condense: 1 } }), { ref: 'c', part: 'kinetochore' }), [-1.2, 0, 0]);
  close(anchorPoint(decl, state({ p: { condense: 1 } }), { ref: 'c', part: 'top' }), [0, 6, 0]);
  close(anchorPoint(decl, state({ pos: [5, 0, 0], p: { condense: 1 } }), { ref: 'c', part: 'centromere' }), [5, 0, 0]);
});

test('tip parts use resolved ends', () => {
  const decl = { type: 'fibre', params: {} };
  close(anchorPoint(decl, state({ ends: { from: [0, 0, 0], to: [10, 0, 0], tip: [5, 0, 0] } }), { ref: 'f', part: 'tip' }), [5, 0, 0]);
});

test('actor extents', () => {
  assert.equal(actorExtent('helix', { turns: 10, scale: 1 }), 17);
  assert.equal(actorExtent('chromosome', { length: 12, scale: 1 }), 6);
  assert.equal(actorExtent('cell', { radius: 40, scale: 2 }), 80);
  assert.equal(actorExtent('ring', { radius: 3, scale: 1 }), 3);
  assert.equal(actorExtent('protein', { size: 6, scale: 1 }), 3);
  assert.equal(actorExtent('point', { scale: 1 }), 1);
});
