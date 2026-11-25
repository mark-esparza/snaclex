import test from 'node:test';
import assert from 'node:assert/strict';
import { parse } from '../parse.js';
import { suggest, VERBS, CAMERA_VERBS, ACTOR_TYPES } from '../schema.js';

const msgs = (r) => r.errors.map((e) => e.message);

test('suggests close words', () => {
  assert.equal(suggest('pul', ['pull', 'push']), 'pull');
  assert.equal(suggest('zzz', ['pull']), null);
  assert.equal(suggest('Helx', ['helix', 'ring']), 'helix');
});

test('unknown verb', () => {
  const r = parse('actor f = fibre from=0,0,0 to=1,0,0\n1s f pul');
  assert.equal(r.scene, null);
  assert.equal(r.errors[0].line, 2);
  assert.match(r.errors[0].message, /did you mean "pull"\?$/);
});

test('unknown header word gets a suggestion', () => {
  const r = parse('titel "x"');
  assert.match(r.errors[0].message, /did you mean "title"\?$/);
});

test('unknown actor type and param key', () => {
  const r = parse('actor d = helx\nactor e = helix trns=4');
  assert.match(msgs(r)[0], /did you mean "helix"\?$/);
  assert.equal(r.errors[1].line, 2);
  assert.match(msgs(r)[1], /did you mean "turns"\?$/);
});

test('verb not for type', () => {
  const r = parse('actor d = helix\n1s d open');
  assert.equal(r.errors[0].message, '"open" does not apply to helix "d"');
});

test('undeclared target with suggestion', () => {
  const r = parse('actor dna = helix\n0-2s dnna unwind from=0.2 to=0.4');
  assert.match(msgs(r)[0], /did you mean "dna"\?/);
});

test('undeclared ref in args with suggestion', () => {
  const r = parse('actor dna = helix\nactor p = protein\n0-1s p move to=dnna@0.4');
  assert.equal(r.errors[0].line, 3);
  assert.match(msgs(r)[0], /did you mean "dna"\?/);
});

test('refs may point at actors declared later', () => {
  const r = parse('actor p = protein at=dna@0.5\nactor dna = helix');
  assert.deepEqual(r.errors, []);
});

test('bad part and @u on a non-helix', () => {
  const r = parse('actor c = chromosome\nactor p = protein\nactor q = protein\n0-1s p move to=c.elbow\n0-1s q move to=p@0.5');
  assert.match(msgs(r)[0], /chromosome "c" has no part "elbow"/);
  assert.match(msgs(r)[1], /only a helix/);
});

test('instant gets default duration', () => {
  const r = parse('actor r = ring\n6s r open');
  assert.deepEqual(r.errors, []);
  assert.equal(r.scene.actions[0].t1, 7);
});

test('caption duration from for= or range', () => {
  const r = parse('2s caption "a" for=4s\n5-6s caption "b" for=9s\n8s caption "c"');
  const [a, b, c] = r.scene.actions;
  assert.equal(a.t1, 6);
  assert.equal(b.t1, 6);
  assert.equal(c.t1, 11);
});

test('positional object maps to the verb argument and bare names become refs', () => {
  const r = parse('actor fibre1 = fibre from=poleL to=chr21a.kinetochore\nactor chr21a = chromosome\nactor poleL = point\n14-20s fibre1 pull chr21a to=poleL\n0s camera zoom chr21a');
  assert.deepEqual(r.errors, []);
  const [pull, zoom] = r.scene.actions;
  assert.deepEqual(pull.args.target, { ref: 'chr21a' });
  assert.deepEqual(pull.args.to, { ref: 'poleL' });
  assert.equal(pull.args.object, undefined);
  assert.deepEqual(zoom.args.to, { ref: 'chr21a' });
  assert.deepEqual(r.scene.actors[0].params.from, { ref: 'poleL' });
});

test('missing required argument', () => {
  const r = parse('actor p = protein\n0-1s p move');
  assert.match(msgs(r)[0], /"move" needs to=/);
});

test('number range checked', () => {
  const r = parse('actor p = protein\n0-1s p fade to=2');
  assert.match(msgs(r)[0], /between 0 and 1/);
});

test('defaults are filled in', () => {
  const r = parse('actor d = helix\nactor p = protein\n1s d unwind\n2s p rotate');
  assert.deepEqual(r.errors, []);
  const d = r.scene.actors[0].params;
  assert.equal(d.turns, 20);
  assert.equal(d.color, '#ff5fa2');
  assert.deepEqual(d.at, [0, 0, 0]);
  assert.equal(d.scale, 1);
  const [unwind, rotate] = r.scene.actions;
  assert.equal(unwind.args.from, 0.45);
  assert.equal(unwind.args.ease, 'smooth');
  assert.equal(rotate.args.axis, 'y');
});

test('duration cap', () => {
  assert.match(msgs(parse('duration 301s'))[0], /300s/);
});

test('bad size', () => {
  const m = msgs(parse('size 800x600'))[0];
  assert.match(m, /1280x720/);
  assert.match(m, /1920x1080/);
});

test('bad style', () => {
  assert.match(msgs(parse('style neon'))[0], /dark/);
});

test('duplicate actor', () => {
  const r = parse('actor a = ring\nactor a = helix');
  assert.equal(r.errors[0].line, 2);
  assert.match(msgs(r)[0], /already declared on line 1/);
});

test('reserved actor names', () => {
  assert.match(msgs(parse('actor camera = ring'))[0], /reserved/);
});

test('protein from must look like a PDB ID', () => {
  assert.match(msgs(parse('actor p = protein from=separase'))[0], /PDB ID/);
  assert.deepEqual(parse('actor p = protein from=7NJ1 chain=A,B').errors, []);
});

test('label needs text', () => {
  assert.match(msgs(parse('actor l = label'))[0], /needs text=/);
});

test('every verb has a default duration and a target list', () => {
  for (const [name, v] of Object.entries(VERBS)) {
    assert.equal(typeof v.defaultDuration, 'number', name);
    assert.ok(Array.isArray(v.appliesTo) && v.appliesTo.length, name);
    for (const t of v.appliesTo) assert.ok(t === '*' || ACTOR_TYPES[t], `${name} → ${t}`);
  }
  for (const [name, v] of Object.entries(CAMERA_VERBS)) {
    assert.equal(typeof v.defaultDuration, 'number', name);
  }
});
