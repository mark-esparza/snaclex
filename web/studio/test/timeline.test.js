import test from 'node:test';
import assert from 'node:assert/strict';
import { parse } from '../parse.js';
import { compile, stateAt, ease } from '../timeline.js';

function build(text) {
  const p = parse(text);
  assert.deepEqual(p.errors, [], 'script should parse');
  const c = compile(p.scene);
  return c;
}
function tl(text) {
  const c = build(text);
  assert.deepEqual(c.errors, [], 'script should compile');
  return c.timeline;
}
const near = (a, b, eps = 1e-6) => assert.ok(Math.abs(a - b) < eps, `${a} ≉ ${b}`);
const nearV = (a, b, eps = 1e-6) => a.forEach((v, i) => near(v, b[i], eps));

test('ease endpoints', () => {
  for (const e of ['smooth', 'linear', 'in', 'out']) {
    assert.equal(ease(e, 0), 0);
    assert.equal(ease(e, 1), 1);
  }
  assert.equal(ease('smooth', 0.5), 0.5);
  assert.equal(ease('linear', 0.25), 0.25);
});

test('move interpolates position', () => {
  const t = tl('actor p = protein at=0,0,0\n0-2s p move to=10,0,0');
  near(stateAt(t, 1).actors.p.pos[0], 5);
  near(stateAt(t, 3).actors.p.pos[0], 10);
  near(stateAt(t, 0).actors.p.pos[0], 0);
});

test('enter hides before start', () => {
  const t = tl('actor p = protein\n2-4s p enter from=left to=0,0,0');
  assert.equal(stateAt(t, 1).actors.p.visible, false);
  const s = stateAt(t, 4).actors.p;
  assert.deepEqual(s.pos, [0, 0, 0]);
  assert.equal(s.opacity, 1);
  assert.equal(s.visible, true);
  assert.ok(stateAt(t, 2.1).actors.p.pos[0] < -40, 'starts off to the left');
});

test('leave ends invisible', () => {
  const t = tl('actor p = protein\n1-2s p leave to=right');
  assert.equal(stateAt(t, 0.5).actors.p.visible, true);
  const s = stateAt(t, 2.5).actors.p;
  assert.equal(s.visible, false);
  assert.equal(s.opacity, 0);
  near(s.pos[0], 60);
});

test('hidden actors and show', () => {
  const t = tl('actor p = protein hidden=yes\n3-4s p show');
  assert.equal(stateAt(t, 1).actors.p.visible, false);
  assert.equal(stateAt(t, 4).actors.p.visible, true);
  assert.equal(stateAt(t, 4).actors.p.opacity, 1);
});

test('ref resolves to helix point', () => {
  const t = tl('actor dna = helix turns=10\nactor p = protein\n0-1s p move to=dna@1');
  nearV(stateAt(t, 1).actors.p.pos, [17, 0, 0]);
});

test('at= may reference another actor', () => {
  const t = tl('actor p = protein at=dna@0\nactor dna = helix turns=10 at=0,5,0');
  nearV(stateAt(t, 0).actors.p.pos, [-17, 5, 0]);
});

test('attach follows a moving target', () => {
  const t = tl('actor p = protein at=2,0,0\nactor q = protein at=0,0,0\n0-0.5s p attach to=q\n1-2s q move to=0,10,0');
  nearV(stateAt(t, 2).actors.p.pos, [2, 10, 0]);
  nearV(stateAt(t, 1.5).actors.p.pos, [2, 5, 0]);
});

test('detach freezes the position', () => {
  const t = tl('actor p = protein at=2,0,0\nactor q = protein\n0s p attach to=q\n1-2s q move to=0,10,0\n3s p detach\n4-5s q move to=0,-10,0');
  nearV(stateAt(t, 6).actors.p.pos, [2, 10, 0]);
});

test('move after attach starts from the attached position', () => {
  const t = tl('actor p = protein at=2,0,0\nactor q = protein\n0s p attach to=q\n1-2s q move to=0,10,0\n3-4s p move to=0,0,0');
  nearV(stateAt(t, 3).actors.p.pos, [2, 10, 0]);
  nearV(stateAt(t, 4).actors.p.pos, [0, 0, 0]);
});

test('ring around follows its target', () => {
  const t = tl('actor dna = helix turns=10\nactor r = ring around=dna@1\n0-1s dna move to=0,4,0');
  nearV(stateAt(t, 1).actors.r.pos, [17, 4, 0]);
  nearV(stateAt(t, 1).actors.r.ends.axis, [1, 0, 0]);
});

test('fibre ends resolve every frame; grow from zero when grown', () => {
  const t = tl('actor pole = point at=-30,0,0\nactor c = chromosome\nactor f = fibre from=pole to=c.kinetochore\n0-2s f grow\n3-4s c move to=10,0,0');
  const s0 = stateAt(t, 0).actors.f;
  assert.equal(s0.p.grow, 0);
  nearV(s0.ends.from, [-30, 0, 0]);
  nearV(s0.ends.to, [-1.2, 0, 0]);
  const s4 = stateAt(t, 4).actors.f;
  nearV(s4.ends.to, [8.8, 0, 0]);
  nearV(s4.ends.tip, [8.8, 0, 0]);
});

test('pull moves the target most of the way', () => {
  const t = tl('actor pole = point at=-30,0,0\nactor c = chromosome\nactor f = fibre from=pole to=c.kinetochore\n0-3s f pull c to=pole');
  nearV(stateAt(t, 3).actors.c.pos, [-21, 0, 0]);
});

test('unwind opens a bubble from its centre, then slides', () => {
  const t = tl('actor d = helix\n0-2s d unwind from=0.4 to=0.5\n2-4s d unwind from=0.5 to=0.6\n4-5s d rewind');
  const a = stateAt(t, 0).actors.d.p;
  near(a.unwind0, 0.45); near(a.unwind1, 0.45);
  const b = stateAt(t, 2).actors.d.p;
  near(b.unwind0, 0.4); near(b.unwind1, 0.5);
  const c = stateAt(t, 4).actors.d.p;
  near(c.unwind0, 0.5); near(c.unwind1, 0.6);
  const d = stateAt(t, 5).actors.d.p;
  near(d.unwind0, d.unwind1);
});

test('state props for ring, chromosome, cell and cut', () => {
  const t = tl('actor r = ring\nactor c = chromosome\nactor k = cell\nactor d = helix\n0-1s r open\n1-2s r cleave\n0-2s c separate gap=8\n0-2s k pinch\n2-4s k split\n0-1s d cut at=0.3');
  const s = stateAt(t, 5).actors;
  assert.equal(s.r.p.open, 1);
  assert.equal(s.r.p.cleave, 1);
  assert.equal(s.c.p.separate, 1);
  assert.equal(s.c.p.gap, 8);
  assert.equal(s.k.p.pinch, 1);
  assert.equal(s.k.p.split, 1);
  assert.equal(s.d.p.cut, 1);
  assert.equal(s.d.p.cutAt, 0.3);
});

test('colour, fade, glow, scale and rotate tween', () => {
  const t = tl('actor p = protein color=#000000\n0-1s p color to=#ffffff ease=linear\n0-1s p fade to=0.5 ease=linear\n0-1s p glow to=1\n0-1s p scale to=3\n0-1s p rotate by=90 axis=z');
  const mid = stateAt(t, 0.5).actors.p;
  assert.equal(mid.color, '#808080');
  near(mid.opacity, 0.75);
  const end = stateAt(t, 1).actors.p;
  assert.equal(end.glow, 1);
  assert.equal(end.scale, 3);
  nearV(end.rot, [0, 0, 90]);
});

test('pair places a homolog beside its partner', () => {
  const t = tl('actor a = chromosome at=0,0,0\nactor b = chromosome at=20,0,0\n0-2s b pair with=a gap=3');
  nearV(stateAt(t, 2).actors.b.pos, [3, 0, 0]);
});

test('unknown ref', () => {
  const p = parse('actor dna = helix\nactor p = protein\n0-1s p move to=dnna@0.4');
  const all = p.scene ? compile(p.scene).errors : p.errors;
  assert.match(all[0].message, /did you mean "dna"\?/);
});

test('attach cycle', () => {
  const c = build('actor a = protein\nactor b = protein\n0s a attach to=b\n1s b attach to=a');
  assert.match(c.errors[0].message, /^attach cycle: /);
  assert.equal(c.timeline, null);
});

test('overlap warning', () => {
  const c = build('actor p = protein\n0-4s p move to=10,0,0\n2-3s p move to=0,5,0');
  assert.equal(c.warnings.length, 1);
  assert.match(c.warnings[0].message, /overrides line \d+ \(pos of "p"\)/);
  assert.equal(c.warnings[0].line, 3);
});

test('duration overflow', () => {
  const c = build('duration 5s\nactor p = protein\n4-6s p move to=1,0,0');
  assert.equal(c.errors[0].message, 'action ends at 6s but duration is 5s');
});

test('duration defaults to last action + 1s', () => {
  assert.equal(tl('credits off\nactor p = protein\n2-6s p move to=1,0,0').duration, 7);
  assert.equal(tl('actor p = protein\n2-6s p move to=1,0,0').duration, 7);
});

test('credits text', () => {
  const t = tl('actor sep = protein from=7NJ1 chain=A label="separase"\nactor sec = protein from=7NJ1 chain=B\n0-2s sep move to=1,0,0');
  assert.ok(t.creditsText.startsWith('Structures: 7NJ1 chain A (separase), 7NJ1 chain B'));
  assert.ok(t.creditsText.endsWith(' · Rendered with SnaCleX Studio'));
  assert.equal(t.duration, 7);
  assert.equal(stateAt(t, t.duration - 1).creditsVisible, true);
  assert.equal(stateAt(t, 1).creditsVisible, false);
  assert.equal(tl('credits off\nactor sep = protein from=7NJ1').creditsText, null);
});

test('caption visibility', () => {
  const t = tl('actor p = protein\n2s caption "Hi" for=3s\n0-6s p move to=1,0,0');
  assert.deepEqual(stateAt(t, 4).captions, ['Hi']);
  assert.deepEqual(stateAt(t, 5.5).captions, []);
});

test('default camera frames the scene', () => {
  const t = tl('actor dna = helix turns=10\nactor p = protein at=0,10,0');
  const cam = stateAt(t, 0).camera;
  assert.equal(cam.fov, 35);
  nearV(cam.target, [0, 5.3, 0]); // centre of the bounding box: helix ±2.4, protein 10 ± 3
  const d = Math.hypot(...cam.pos.map((v, i) => v - cam.target[i]));
  assert.ok(d > 40 && d < 200, `distance ${d}`);
  assert.ok(cam.pos[2] > cam.target[2], 'looks from +z');
});

test('camera zoom, orbit, follow and cut', () => {
  const t = tl('actor a = protein at=0,0,0 size=4\nactor b = protein at=50,0,0\n0-2s camera zoom to=a distance=10\n2-4s camera orbit around=a degrees=180 ease=linear\n4-6s camera follow b\n4-6s b move to=50,20,0\n7s camera cut to=a distance=5');
  const z = stateAt(t, 2).camera;
  nearV(z.target, [0, 0, 0]);
  near(Math.hypot(...z.pos), 10);
  const o = stateAt(t, 4).camera;
  near(o.pos[2], -z.pos[2], 1e-6);
  const f = stateAt(t, 6).camera;
  nearV(f.target, [50, 20, 0]);
  const c = stateAt(t, 7).camera;
  nearV(c.target, [0, 0, 0]);
  near(Math.hypot(...c.pos), 5);
});

test('stateAt is order-independent', () => {
  const t = tl('actor dna = helix\nactor p = protein\nactor q = protein at=5,5,0\n0-2s p enter to=dna@0.3\n1-3s dna unwind\n2s p attach to=q\n3-5s q move to=0,10,0\n4s caption "x"\n5-7s camera zoom p');
  const first = JSON.stringify(stateAt(t, 3));
  stateAt(t, 9); stateAt(t, 0); stateAt(t, 5);
  assert.equal(JSON.stringify(stateAt(t, 3)), first);
});

test('stateAt does not mutate the timeline', () => {
  const t = tl('actor p = protein\n0-2s p move to=10,0,0');
  const before = JSON.stringify(t);
  const s = stateAt(t, 1);
  s.actors.p.pos[0] = 999;
  stateAt(t, 2);
  assert.equal(JSON.stringify(t), before);
});

test('wide framing fits a long helix to the 16:9 frame width', () => {
  const t = tl('actor dna = helix turns=20');
  const cam = stateAt(t, 0).camera;
  const d = Math.hypot(...cam.pos.map((v, i) => v - cam.target[i]));
  assert.ok(d > 55 && d < 80, `distance ${d}`);
});
