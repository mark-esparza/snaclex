import test from 'node:test';
import assert from 'node:assert/strict';
import { parse, parseSyntax } from '../parse.js';

const SKETCH = `title "Nondisjunction in Meiosis I"
size 1920x1080
duration 45s
style dark
note "Stylised; not to scale"

# --- actors ---
actor dna      = helix turns=40 color=pink
actor sep      = protein from=7NJ1 chain=A color=gold
actor sec      = protein from=7NJ1 chain=B color=orange
actor ring     = ring color=teal around=dna at=0.5
actor chr21a   = chromosome color=violet label="Chr 21 (maternal)"
actor poleL    = point at=-30,0,0
actor fibre1   = fibre from=poleL to=chr21a.kinetochore

# --- timeline ---
0s      camera wide
0-4s    sep enter from=left to=dna@0.45
4-6s    sec leave to=right
6s      caption "Securin released: separase is now active" for=3s
6-8s    ring open
8-12s   camera zoom to=ring
14-20s  fibre1 pull chr21a to=poleL
`;

test('parses header', () => {
  const r = parseSyntax('title "X"\nsize 1280x720\nduration 45s\nstyle dark\nnote "n"\ncredits off');
  assert.deepEqual(r.errors, []);
  assert.equal(r.scene.title, 'X');
  assert.deepEqual(r.scene.size, [1280, 720]);
  assert.equal(r.scene.duration, 45);
  assert.equal(r.scene.style, 'dark');
  assert.equal(r.scene.note, 'n');
  assert.equal(r.scene.credits, false);
});

test('header defaults', () => {
  const r = parseSyntax('');
  assert.deepEqual(r.scene.size, [1920, 1080]);
  assert.equal(r.scene.duration, null);
  assert.equal(r.scene.style, 'dark');
  assert.equal(r.scene.credits, true);
});

test('parses actor with typed values', () => {
  const r = parseSyntax('actor sep = protein from=7NJ1 chain=A color=gold at=1,2,3');
  assert.deepEqual(r.errors, []);
  const a = r.scene.actors[0];
  assert.equal(a.name, 'sep');
  assert.equal(a.type, 'protein');
  assert.equal(a.line, 1);
  assert.deepEqual(a.params, { from: '7NJ1', chain: 'A', color: '#ffc94d', at: [1, 2, 3] });
});

test('parses refs, dirs, hex colours, booleans and quoted strings', () => {
  const r = parseSyntax('actor x = protein a=dna@0.45 b=chr21a.kinetochore c=left d=#AbCdEf e=yes f="two words" g=-1.5,0,2e1');
  assert.deepEqual(r.errors, []);
  assert.deepEqual(r.scene.actors[0].params, {
    a: { ref: 'dna', u: 0.45 },
    b: { ref: 'chr21a', part: 'kinetochore' },
    c: { dir: 'left' },
    d: '#abcdef',
    e: true,
    f: 'two words',
    g: [-1.5, 0, 20],
  });
});

test('parses ranges and instants', () => {
  const r = parseSyntax('0-4s sep enter from=left\n6s ring open\n2.5-3.25s x move\n0s-4s y move\n7 z show');
  assert.deepEqual(r.errors, []);
  const [a, b, c, d, e] = r.scene.actions;
  assert.equal(a.t0, 0); assert.equal(a.t1, 4); assert.equal(a.verb, 'enter'); assert.equal(a.target, 'sep');
  assert.deepEqual(a.args, { from: { dir: 'left' } });
  assert.equal(b.t0, 6); assert.equal(b.t1, null); assert.equal(b.verb, 'open');
  assert.equal(c.t0, 2.5); assert.equal(c.t1, 3.25);
  assert.equal(d.t0, 0); assert.equal(d.t1, 4);
  assert.equal(e.t0, 7); assert.equal(e.t1, null);
});

test('parses positional object', () => {
  const r = parseSyntax('14-20s fibre1 pull chr21a to=poleL');
  assert.deepEqual(r.errors, []);
  const a = r.scene.actions[0];
  assert.equal(a.verb, 'pull');
  assert.deepEqual(a.args.object, { ref: 'chr21a' });
  assert.equal(a.args.to, 'poleL');
});

test('camera lines', () => {
  const r = parseSyntax('0s camera wide\n1-3s camera zoom ring');
  assert.deepEqual(r.errors, []);
  assert.equal(r.scene.actions[0].target, 'camera');
  assert.equal(r.scene.actions[0].verb, 'wide');
  assert.deepEqual(r.scene.actions[1].args.object, { ref: 'ring' });
});

test('caption with text', () => {
  const r = parseSyntax('6s caption "Securin released" for=3s');
  assert.deepEqual(r.errors, []);
  const a = r.scene.actions[0];
  assert.equal(a.target, 'caption');
  assert.equal(a.verb, 'caption');
  assert.equal(a.text, 'Securin released');
  assert.equal(a.args.for, 3);
});

test('comments and quotes', () => {
  const r = parseSyntax('title "A # B" # c\n# whole line\n  \nactor d = helix color=#ff0000 # trailing');
  assert.deepEqual(r.errors, []);
  assert.equal(r.scene.title, 'A # B');
  assert.equal(r.scene.actors[0].params.color, '#ff0000');
});

test('escaped quote inside a string', () => {
  const r = parseSyntax('title "say \\"hi\\""');
  assert.equal(r.scene.title, 'say "hi"');
});

test('reversed range errors with line', () => {
  const r = parseSyntax('title "a"\n\n5-2s x move');
  assert.equal(r.scene, null);
  assert.deepEqual(r.errors[0], { line: 3, message: 'end time is before start time' });
});

test('unterminated string', () => {
  const r = parseSyntax('title "oops');
  assert.equal(r.errors.length, 1);
  assert.match(r.errors[0].message, /missing closing "/);
});

test('malformed lines report their line number', () => {
  const r = parseSyntax('actor = helix\nactor d helix\nsize big\n0s\n0s caption nope\n1s a move 3,4');
  const lines = r.errors.map((e) => e.line);
  assert.deepEqual(lines, [1, 2, 3, 4, 5, 6]);
});

test('unknown first word is an error', () => {
  const r = parseSyntax('wobble 3');
  assert.equal(r.errors[0].line, 1);
  assert.match(r.errors[0].message, /"wobble"/);
});

test('the spec sketch parses without syntax errors', () => {
  const r = parseSyntax(SKETCH);
  assert.deepEqual(r.errors, []);
  assert.equal(r.scene.actors.length, 7);
  assert.equal(r.scene.actions.length, 7);
});

test('never throws on truncated input', () => {
  let seed = 7;
  const rand = () => ((seed = (seed * 1103515245 + 12345) % 2147483648) / 2147483648);
  for (let i = 0; i < 200; i++) {
    const cut = SKETCH.slice(0, Math.floor(rand() * SKETCH.length));
    const r = parse(cut);
    assert.equal(typeof r, 'object');
    assert.ok(Array.isArray(r.errors));
  }
});
