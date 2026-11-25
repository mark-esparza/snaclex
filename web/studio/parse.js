// SnaCleX Studio — scene-script parser.
//
// The script is line based. Each non-blank line is a header (`title "…"`), an
// actor declaration (`actor dna = helix turns=40`) or a timeline action
// (`0-4s sep enter from=left`). parseSyntax() types each value by how it is
// written; parse() also checks it against the schema (actor types, verbs,
// parameter kinds). Neither ever throws: problems are collected as
// {line, message} diagnostics so a half-typed script still gives useful
// feedback.
//
// Pure module: no DOM, no WebGL. Node's test runner imports it directly.

export const COLOR_NAMES = {
  pink: '#ff5fa2', red: '#ff4d4d', orange: '#ff9a3c', gold: '#ffc94d',
  yellow: '#fff06b', green: '#5fe07a', teal: '#3fd6c6', blue: '#4da3ff',
  violet: '#a77bff', white: '#f2f2f2', grey: '#9aa0a6',
};

export const DIRS = ['left', 'right', 'top', 'bottom', 'front', 'back'];
export const HEADER_KEYS = ['title', 'size', 'duration', 'style', 'note', 'credits'];

const NAME_RE = /^[A-Za-z][A-Za-z0-9_-]*$/;
const NUM = '-?(?:\\d+\\.?\\d*|\\.\\d+)(?:e-?\\d+)?';
const NUM_RE = new RegExp(`^${NUM}s?$`, 'i');
const VEC_RE = new RegExp(`^(${NUM}),(${NUM}),(${NUM})$`, 'i');
const TIME = '(\\d+(?:\\.\\d+)?)s?';
const TIMESPEC_RE = new RegExp(`^${TIME}(?:-${TIME})?$`);
const REF_U_RE = new RegExp(`^([A-Za-z][A-Za-z0-9_-]*)@(${NUM})$`);
const REF_PART_RE = /^([A-Za-z][A-Za-z0-9_-]*)\.([A-Za-z][A-Za-z0-9_-]*)$/;

class LineError extends Error {}

// Split a line into tokens. A token is a run of non-space characters in which
// "…" segments may contain spaces. `#` starts a comment only at the start of a
// token, so `color=#ff5fa2` stays a value. Returns [{text, quoted}] where
// `text` has quotes removed and `quoted` says whether any part was quoted.
export function tokenize(line) {
  const tokens = [];
  let i = 0;
  const n = line.length;
  while (i < n) {
    while (i < n && /\s/.test(line[i])) i++;
    if (i >= n) break;
    if (line[i] === '#') break;
    let text = '';
    let quoted = false;
    while (i < n && !/\s/.test(line[i])) {
      if (line[i] === '"') {
        quoted = true;
        i++;
        let closed = false;
        while (i < n) {
          if (line[i] === '\\' && i + 1 < n) { text += line[i + 1]; i += 2; continue; }
          if (line[i] === '"') { closed = true; i++; break; }
          text += line[i++];
        }
        if (!closed) throw new LineError('missing closing " on this line');
      } else {
        text += line[i++];
      }
    }
    tokens.push({ text, quoted });
  }
  return tokens;
}

// Type a bare (unquoted) value by its spelling.
export function typeValue(text) {
  const low = text.toLowerCase();
  if (low === 'yes' || low === 'on') return true;
  if (low === 'no' || low === 'off') return false;
  if (NUM_RE.test(text)) return parseFloat(text);
  if (/^#[0-9a-f]{6}$/i.test(text)) return text.toLowerCase();
  if (Object.hasOwn(COLOR_NAMES, low)) return COLOR_NAMES[low];
  let m = VEC_RE.exec(text);
  if (m) return [parseFloat(m[1]), parseFloat(m[2]), parseFloat(m[3])];
  m = REF_U_RE.exec(text);
  if (m) return { ref: m[1], u: parseFloat(m[2]) };
  m = REF_PART_RE.exec(text);
  if (m) return { ref: m[1], part: m[2] };
  if (DIRS.includes(low)) return { dir: low };
  return text;
}

function parseRefToken(tok) {
  if (tok.quoted) return null;
  let m = REF_U_RE.exec(tok.text);
  if (m) return { ref: m[1], u: parseFloat(m[2]) };
  m = REF_PART_RE.exec(tok.text);
  if (m) return { ref: m[1], part: m[2] };
  if (NAME_RE.test(tok.text)) return { ref: tok.text };
  return null;
}

function parseKeyValues(tokens, start, into) {
  for (let i = start; i < tokens.length; i++) {
    const tok = tokens[i];
    const eq = tok.text.indexOf('=');
    if (eq <= 0 || (tok.quoted && !/^[A-Za-z][A-Za-z0-9_-]*=/.test(tok.text))) {
      throw new LineError(`expected key=value but found "${tok.text}"`);
    }
    const key = tok.text.slice(0, eq);
    if (!NAME_RE.test(key)) throw new LineError(`"${key}" is not a valid key`);
    const raw = tok.text.slice(eq + 1);
    if (!tok.quoted && raw === '') throw new LineError(`"${key}=" needs a value`);
    into[key] = tok.quoted ? raw : typeValue(raw);
  }
  return into;
}

function parseHeader(key, tokens, scene) {
  const arg = tokens[1];
  const need = (what) => { throw new LineError(`${key} needs ${what}`); };
  if (tokens.length > 2) throw new LineError(`unexpected "${tokens[2].text}" after ${key}`);
  switch (key) {
    case 'title':
    case 'note':
      if (!arg || !arg.quoted) need('a quoted text, like ' + key + ' "…"');
      scene[key] = arg.text;
      return;
    case 'size': {
      const m = arg && !arg.quoted && /^(\d+)x(\d+)$/.exec(arg.text);
      if (!m) need('a size like 1920x1080');
      scene.size = [parseInt(m[1], 10), parseInt(m[2], 10)];
      return;
    }
    case 'duration': {
      const m = arg && !arg.quoted && /^(\d+(?:\.\d+)?)s?$/.exec(arg.text);
      if (!m) need('a time like 45s');
      scene.duration = parseFloat(m[1]);
      return;
    }
    case 'style':
      if (!arg || arg.quoted) need('a style name, like dark');
      scene.style = arg.text.toLowerCase();
      return;
    case 'credits': {
      const v = arg && !arg.quoted ? arg.text.toLowerCase() : '';
      if (v !== 'on' && v !== 'off') need('on or off');
      scene.credits = v === 'on';
      return;
    }
  }
}

function parseActor(tokens, line) {
  const name = tokens[1];
  if (!name || name.quoted || !NAME_RE.test(name.text)) {
    throw new LineError('actor needs a name, like: actor dna = helix');
  }
  if (!tokens[2] || tokens[2].text !== '=' || tokens[2].quoted) {
    throw new LineError(`expected "=" after the actor name, like: actor ${name.text} = helix`);
  }
  const type = tokens[3];
  if (!type || type.quoted || !NAME_RE.test(type.text)) {
    throw new LineError(`actor "${name.text}" needs a type after "="`);
  }
  return { name: name.text, type: type.text.toLowerCase(), params: parseKeyValues(tokens, 4, {}), line };
}

function parseAction(tokens, line, m) {
  const t0 = parseFloat(m[1]);
  const t1 = m[2] !== undefined ? parseFloat(m[2]) : null;
  if (t1 !== null && t1 < t0) throw new LineError('end time is before start time');
  const target = tokens[1];
  if (!target || target.quoted || !NAME_RE.test(target.text)) {
    throw new LineError('after the time, name what acts — an actor, camera or caption');
  }
  if (target.text === 'caption') {
    const text = tokens[2];
    if (!text || !text.quoted) throw new LineError('caption needs quoted text, like: 6s caption "Hello"');
    return { t0, t1, target: 'caption', verb: 'caption', args: parseKeyValues(tokens, 3, {}), text: text.text, line };
  }
  const verb = tokens[2];
  if (!verb || verb.quoted || !NAME_RE.test(verb.text)) {
    throw new LineError(`what should "${target.text}" do? Add a verb, like: move`);
  }
  const args = {};
  let next = 3;
  const obj = tokens[3];
  if (obj && !obj.text.includes('=')) {
    const ref = parseRefToken(obj);
    if (!ref) throw new LineError(`unexpected "${obj.text}"`);
    args.object = ref;
    next = 4;
  }
  parseKeyValues(tokens, next, args);
  return { t0, t1, target: target.text, verb: verb.text.toLowerCase(), args, text: null, line };
}

// Syntax only: build the Scene shape with values typed by spelling.
export function parseSyntax(text) {
  const scene = {
    title: '', size: [1920, 1080], duration: null, style: 'dark',
    note: null, credits: true, actors: [], actions: [],
  };
  const errors = [];
  const lines = String(text ?? '').split(/\r?\n/);
  lines.forEach((raw, idx) => {
    const line = idx + 1;
    try {
      const tokens = tokenize(raw);
      if (tokens.length === 0) return;
      const first = tokens[0];
      const word = first.quoted ? '' : first.text;
      if (HEADER_KEYS.includes(word)) return parseHeader(word, tokens, scene);
      if (word === 'actor') return void scene.actors.push(parseActor(tokens, line));
      const m = TIMESPEC_RE.exec(word);
      if (m) return void scene.actions.push(parseAction(tokens, line, m));
      throw new LineError(unknownStart(first.text));
    } catch (e) {
      if (!(e instanceof LineError)) throw e;
      errors.push({ line, message: e.message });
    }
  });
  return { scene: errors.length ? null : scene, errors, warnings: [] };
}

// Hook so the schema module can add "did you mean" to unknown first words.
let unknownStartHook = null;
export function setUnknownStartHook(fn) { unknownStartHook = fn; }
function unknownStart(word) {
  const base = `"${word}" is not a header, "actor", or a time`;
  const hint = unknownStartHook ? unknownStartHook(word) : null;
  return hint ? `${base} — did you mean "${hint}"?` : base;
}

// Full parse: syntax, then schema validation (added by the schema module).
let validator = null;
export function setValidator(fn) { validator = fn; }

export function parse(text) {
  const r = parseSyntax(text);
  if (!r.scene || !validator) return r;
  const v = validator(r.scene);
  return {
    scene: v.errors.length ? null : r.scene,
    errors: v.errors,
    warnings: v.warnings,
  };
}
