// SnaCleX Studio — the scene-script vocabulary and its validation.
//
// Everything a script may say lives in the tables below: header words, actor
// types (with their parameters, named parts and animated state), and verbs
// (with their arguments and default durations). validate() checks a parsed
// Scene against them, coerces bare words into references where a reference is
// expected, fills in defaults, and gives instant actions their default length.
//
// Pure module: no DOM, no WebGL.

export const COLOR_NAMES = {
  pink: '#ff5fa2', red: '#ff4d4d', orange: '#ff9a3c', gold: '#ffc94d',
  yellow: '#fff06b', green: '#5fe07a', teal: '#3fd6c6', blue: '#4da3ff',
  violet: '#a77bff', white: '#f2f2f2', grey: '#9aa0a6',
};
export const DIRS = ['left', 'right', 'top', 'bottom', 'front', 'back'];
export const HEADER_KEYS = ['title', 'size', 'duration', 'style', 'note', 'credits'];
export const SIZES = [[1280, 720], [1920, 1080]];
export const STYLES = ['dark', 'light'];
export const MAX_DURATION = 300;
export const EASINGS = ['smooth', 'linear', 'in', 'out'];
export const RESERVED_NAMES = ['camera', 'caption'];

const unit = { kind: 'number', min: 0, max: 1 };

// Parameters every actor accepts.
const COMMON_PARAMS = {
  at: { kind: 'refOrVec', default: [0, 0, 0] },
  color: { kind: 'color' },
  scale: { kind: 'number', min: 0.01, max: 100, default: 1 },
  rot: { kind: 'vec3', default: [0, 0, 0] },
  label: { kind: 'string', default: null },
  hidden: { kind: 'bool', default: false },
};

export const ACTOR_TYPES = {
  helix: {
    color: '#ff5fa2',
    params: { turns: { kind: 'number', min: 1, max: 200, default: 20 }, curve: { ...unit, default: 0 } },
    parts: ['start', 'end'],
    stateProps: { unwind0: 0, unwind1: 0, cut: 0, cutAt: 0.5, gap: 4 },
  },
  strand: {
    color: '#4da3ff',
    params: { length: { kind: 'number', min: 1, max: 500, default: 10 }, along: { kind: 'ref' } },
    parts: ['tip'],
    stateProps: { grow: 1, peel: 0 },
  },
  protein: {
    color: '#ffc94d',
    params: {
      from: { kind: 'pdbid' },
      chain: { kind: 'string', default: 'A' },
      size: { kind: 'number', min: 0.5, max: 100, default: 6 },
    },
    parts: [],
    stateProps: {},
  },
  ring: {
    color: '#3fd6c6',
    params: { around: { kind: 'ref' }, radius: { kind: 'number', min: 0.2, max: 100, default: 3 } },
    parts: [],
    stateProps: { open: 0, cleave: 0 },
  },
  chromosome: {
    color: '#a77bff',
    params: {
      sisters: { kind: 'number', min: 1, max: 2, default: 2 },
      length: { kind: 'number', min: 2, max: 200, default: 12 },
    },
    parts: ['centromere', 'kinetochore', 'top', 'bottom'],
    stateProps: { separate: 0, condense: 1, gap: 6 },
  },
  fibre: {
    color: '#9fe3ff',
    params: { from: { kind: 'refOrVec', required: true }, to: { kind: 'refOrVec', required: true } },
    parts: ['tip'],
    stateProps: { grow: 1 },
  },
  cell: {
    color: '#9aa0a6',
    params: { radius: { kind: 'number', min: 1, max: 1000, default: 40 }, nucleus: { kind: 'bool', default: true } },
    parts: [],
    stateProps: { pinch: 0, split: 0, envelope: 1, gap: 50 },
  },
  point: {
    color: '#f2f2f2',
    params: { show: { kind: 'bool', default: false } },
    parts: [],
    stateProps: {},
  },
  label: {
    color: '#f2f2f2',
    params: { text: { kind: 'string', required: true }, on: { kind: 'ref' } },
    parts: [],
    stateProps: {},
  },
  arrow: {
    color: '#f2f2f2',
    params: { from: { kind: 'refOrVec', required: true }, to: { kind: 'refOrVec', required: true } },
    parts: [],
    stateProps: {},
  },
};

const ANY = ['*'];
export const VERBS = {
  enter: { appliesTo: ANY, args: { from: { kind: 'refOrDir', default: { dir: 'left' } }, to: { kind: 'refOrVec' } }, defaultDuration: 1.5 },
  leave: { appliesTo: ANY, args: { to: { kind: 'refOrDir', default: { dir: 'right' } } }, defaultDuration: 1.5 },
  move: { appliesTo: ANY, args: { to: { kind: 'refOrVec', required: true } }, defaultDuration: 1.5 },
  attach: { appliesTo: ANY, args: { to: { kind: 'ref', required: true } }, object: 'to', defaultDuration: 0.5 },
  detach: { appliesTo: ANY, args: {}, defaultDuration: 0 },
  rotate: { appliesTo: ANY, args: { by: { kind: 'number', default: 360 }, axis: { kind: 'enum', values: ['x', 'y', 'z'], default: 'y' } }, defaultDuration: 1.5 },
  scale: { appliesTo: ANY, args: { to: { kind: 'number', min: 0, max: 100, required: true } }, defaultDuration: 1 },
  fade: { appliesTo: ANY, args: { to: { ...unit, required: true } }, defaultDuration: 1 },
  color: { appliesTo: ANY, args: { to: { kind: 'color', required: true } }, defaultDuration: 1 },
  glow: { appliesTo: ANY, args: { to: { ...unit, required: true } }, defaultDuration: 1 },
  show: { appliesTo: ANY, args: {}, defaultDuration: 0.5 },
  hide: { appliesTo: ANY, args: {}, defaultDuration: 0.5 },
  unwind: { appliesTo: ['helix'], args: { from: { ...unit, default: 0.45 }, to: { ...unit, default: 0.55 } }, defaultDuration: 2 },
  rewind: { appliesTo: ['helix'], args: {}, defaultDuration: 2 },
  cut: { appliesTo: ['helix'], args: { at: { ...unit, default: 0.5 }, gap: { kind: 'number', min: 0, max: 500, default: 4 } }, defaultDuration: 1.5 },
  grow: { appliesTo: ['strand', 'fibre'], args: { to: { ...unit, default: 1 } }, defaultDuration: 2 },
  peel: { appliesTo: ['strand'], args: {}, defaultDuration: 1.5 },
  open: { appliesTo: ['ring'], args: {}, defaultDuration: 1 },
  close: { appliesTo: ['ring'], args: {}, defaultDuration: 1 },
  cleave: { appliesTo: ['ring'], args: {}, defaultDuration: 1 },
  separate: { appliesTo: ['chromosome'], args: { gap: { kind: 'number', min: 0, max: 500, default: 6 } }, defaultDuration: 2 },
  condense: { appliesTo: ['chromosome'], args: { to: { ...unit, default: 1 } }, defaultDuration: 2 },
  pair: { appliesTo: ['chromosome'], args: { with: { kind: 'ref', required: true }, gap: { kind: 'number', min: 0, max: 500, default: 3 } }, object: 'with', defaultDuration: 2 },
  shrink: { appliesTo: ['fibre'], args: {}, defaultDuration: 2 },
  pull: { appliesTo: ['fibre'], args: { target: { kind: 'ref', required: true }, to: { kind: 'refOrVec', required: true } }, object: 'target', defaultDuration: 3 },
  pinch: { appliesTo: ['cell'], args: { to: { ...unit, default: 1 } }, defaultDuration: 2 },
  split: { appliesTo: ['cell'], args: { gap: { kind: 'number', min: 0, max: 1000, default: 50 } }, defaultDuration: 3 },
};

export const CAMERA_VERBS = {
  wide: { args: {}, defaultDuration: 1.5 },
  zoom: { args: { to: { kind: 'ref', required: true }, distance: { kind: 'number', min: 0.1, max: 10000 } }, object: 'to', defaultDuration: 2 },
  orbit: { args: { around: { kind: 'ref' }, degrees: { kind: 'number', default: 360 } }, object: 'around', defaultDuration: 6 },
  follow: { args: { ref: { kind: 'ref', required: true } }, object: 'ref', defaultDuration: 4 },
  cut: { args: { to: { kind: 'ref', required: true }, distance: { kind: 'number', min: 0.1, max: 10000 } }, object: 'to', defaultDuration: 0 },
};

export const CAPTION_ARGS = { for: { kind: 'number', min: 0.1, max: MAX_DURATION, default: 3 } };
const EASE_ARG = { ease: { kind: 'enum', values: EASINGS, default: 'smooth' } };

// ---------------------------------------------------------------------------

export function levenshtein(a, b) {
  const m = a.length, n = b.length;
  let prev = Array.from({ length: n + 1 }, (_, j) => j);
  for (let i = 1; i <= m; i++) {
    const cur = [i];
    for (let j = 1; j <= n; j++) {
      cur[j] = Math.min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (a[i - 1] === b[j - 1] ? 0 : 1));
    }
    prev = cur;
  }
  return prev[n];
}

// The closest candidate within edit distance 2 (case-insensitive), or null.
export function suggest(word, candidates) {
  const w = String(word).toLowerCase();
  let best = null, bestD = 3;
  for (const c of candidates) {
    const d = levenshtein(w, c.toLowerCase());
    if (d < bestD) { best = c; bestD = d; }
  }
  return best;
}

const didYouMean = (word, candidates) => {
  const s = suggest(word, candidates);
  return s ? ` — did you mean "${s}"?` : '';
};

const isRef = (v) => v && typeof v === 'object' && !Array.isArray(v) && typeof v.ref === 'string';
const isVec = (v) => Array.isArray(v) && v.length === 3 && v.every(Number.isFinite);
const isDir = (v) => v && typeof v === 'object' && typeof v.dir === 'string';
const NAME_RE = /^[A-Za-z][A-Za-z0-9_-]*$/;

const KIND_TEXT = {
  number: 'a number', string: 'text', color: 'a colour (a name like pink, or #rrggbb)',
  vec3: 'a position like 1,2,3', ref: 'the name of an actor', dir: 'a direction',
  refOrVec: 'an actor (like dna@0.5) or a position like 1,2,3',
  refOrDir: 'an actor, a position, or a direction (left, right, top, bottom, front, back)',
  bool: 'yes or no', pdbid: 'a 4-character PDB ID like 7NJ1',
};

// Coerce one value to its spec. Returns {value} or {error}.
function coerce(key, v, spec) {
  const bad = () => ({ error: spec.kind === 'enum'
    ? `"${key}" should be one of ${spec.values.join(', ')}`
    : `"${key}" should be ${KIND_TEXT[spec.kind]}` });
  switch (spec.kind) {
    case 'number': {
      if (typeof v !== 'number' || !Number.isFinite(v)) return bad();
      if ((spec.min !== undefined && v < spec.min) || (spec.max !== undefined && v > spec.max)) {
        return { error: `"${key}" must be between ${spec.min} and ${spec.max}` };
      }
      return { value: v };
    }
    case 'string':
      if (typeof v === 'number') return { value: String(v) };
      if (isRef(v) && v.u === undefined && v.part === undefined) return { value: v.ref };
      if (typeof v !== 'string') return bad();
      return { value: v };
    case 'pdbid': {
      const s = typeof v === 'number' ? String(v) : v;
      if (typeof s !== 'string' || !/^[0-9][A-Za-z0-9]{3}$/.test(s)) return bad();
      return { value: s.toUpperCase() };
    }
    case 'color':
      return typeof v === 'string' && /^#[0-9a-f]{6}$/.test(v) ? { value: v } : bad();
    case 'vec3':
      return isVec(v) ? { value: v } : bad();
    case 'bool':
      return typeof v === 'boolean' ? { value: v } : bad();
    case 'enum': {
      const s = typeof v === 'string' ? v.toLowerCase() : isRef(v) ? v.ref.toLowerCase() : null;
      return spec.values.includes(s) ? { value: s } : bad();
    }
    case 'ref':
      if (isRef(v)) return { value: v };
      if (typeof v === 'string' && NAME_RE.test(v)) return { value: { ref: v } };
      return bad();
    case 'refOrVec':
      if (isVec(v)) return { value: v };
      return coerce(key, v, { kind: 'ref' }).value ? coerce(key, v, { kind: 'ref' }) : bad();
    case 'refOrDir':
      if (isDir(v) || isVec(v)) return { value: v };
      return coerce(key, v, { kind: 'ref' }).value ? coerce(key, v, { kind: 'ref' }) : bad();
    default:
      return bad();
  }
}

// Check and coerce a bag of key=values against a spec table, filling defaults.
function checkBag(bag, specs, what, line, errors) {
  for (const key of Object.keys(bag)) {
    if (!Object.hasOwn(specs, key)) {
      errors.push({ line, message: `${what} has no setting "${key}"${didYouMean(key, Object.keys(specs))}` });
      delete bag[key];
      continue;
    }
    const r = coerce(key, bag[key], specs[key]);
    if (r.error) { errors.push({ line, message: r.error }); delete bag[key]; } else bag[key] = r.value;
  }
  for (const [key, spec] of Object.entries(specs)) {
    if (Object.hasOwn(bag, key)) continue;
    if (spec.required) continue; // reported by the caller with verb-specific wording
    if (spec.default !== undefined) {
      bag[key] = Array.isArray(spec.default) ? [...spec.default]
        : spec.default && typeof spec.default === 'object' ? { ...spec.default } : spec.default;
    }
  }
}

function collectRefs(bag) {
  return Object.values(bag).filter(isRef);
}

export function validate(scene) {
  const errors = [];
  const warnings = [];
  const hl = (key) => (scene.lines && scene.lines[key]) || 1;

  if (scene.duration !== null && scene.duration > MAX_DURATION) {
    errors.push({ line: hl('duration'), message: `maximum duration is ${MAX_DURATION}s` });
  }
  if (!SIZES.some(([w, h]) => w === scene.size[0] && h === scene.size[1])) {
    errors.push({ line: hl('size'), message: `size must be ${SIZES.map(([w, h]) => `${w}x${h}`).join(' or ')}` });
  }
  if (!STYLES.includes(scene.style)) {
    errors.push({ line: hl('style'), message: `style must be ${STYLES.join(' or ')}${didYouMean(scene.style, STYLES)}` });
  }

  // Actors: names, types, params.
  const byName = new Map();
  for (const a of scene.actors) {
    if (RESERVED_NAMES.includes(a.name)) {
      errors.push({ line: a.line, message: `"${a.name}" is reserved; choose another actor name` });
      continue;
    }
    if (byName.has(a.name)) {
      errors.push({ line: a.line, message: `actor "${a.name}" is already declared on line ${byName.get(a.name).line}` });
      continue;
    }
    const type = ACTOR_TYPES[a.type];
    if (!type) {
      errors.push({ line: a.line, message: `unknown actor type "${a.type}"${didYouMean(a.type, Object.keys(ACTOR_TYPES))}` });
      continue;
    }
    byName.set(a.name, a);
    const specs = { ...COMMON_PARAMS, ...type.params };
    checkBag(a.params, specs, `a ${a.type}`, a.line, errors);
    if (a.params.color === undefined) a.params.color = type.color;
    for (const [key, spec] of Object.entries(type.params)) {
      if (spec.required && a.params[key] === undefined) {
        errors.push({ line: a.line, message: `${a.type} "${a.name}" needs ${key}=` });
      }
    }
  }

  const checkRef = (ref, line) => {
    const target = byName.get(ref.ref);
    if (!target) {
      errors.push({ line, message: `no actor named "${ref.ref}"${didYouMean(ref.ref, [...byName.keys()])}` });
      return;
    }
    if (ref.part !== undefined && !ACTOR_TYPES[target.type].parts.includes(ref.part)) {
      const parts = ACTOR_TYPES[target.type].parts;
      errors.push({ line, message: `${target.type} "${target.name}" has no part "${ref.part}"${parts.length ? ` (parts: ${parts.join(', ')})` : ''}` });
    }
    if (ref.u !== undefined) {
      if (target.type !== 'helix') errors.push({ line, message: `"${ref.ref}@${ref.u}": only a helix takes @position` });
      else if (ref.u < 0 || ref.u > 1) errors.push({ line, message: `"${ref.ref}@${ref.u}": position must be between 0 and 1` });
    }
  };
  for (const a of byName.values()) for (const r of collectRefs(a.params)) checkRef(r, a.line);

  // Actions.
  for (const act of scene.actions) {
    const errCount = errors.length;
    let table;
    if (act.target === 'caption') {
      checkBag(act.args, CAPTION_ARGS, 'a caption', act.line, errors);
      if (act.t1 === null) act.t1 = act.t0 + act.args.for;
      continue;
    } else if (act.target === 'camera') {
      table = CAMERA_VERBS[act.verb];
      if (!table) {
        errors.push({ line: act.line, message: `the camera can't "${act.verb}"${didYouMean(act.verb, Object.keys(CAMERA_VERBS))}` });
        continue;
      }
    } else {
      const actor = byName.get(act.target);
      if (!actor) {
        if (!scene.actors.some((a) => a.name === act.target)) {
          errors.push({ line: act.line, message: `no actor named "${act.target}"${didYouMean(act.target, [...byName.keys(), 'camera', 'caption'])}` });
        }
        continue;
      }
      table = VERBS[act.verb];
      const applies = (v) => v.appliesTo.includes('*') || v.appliesTo.includes(actor.type);
      if (!table) {
        const options = Object.keys(VERBS).filter((k) => applies(VERBS[k]));
        errors.push({ line: act.line, message: `unknown verb "${act.verb}"${didYouMean(act.verb, options)}` });
        continue;
      }
      if (!applies(table)) {
        errors.push({ line: act.line, message: `"${act.verb}" does not apply to ${actor.type} "${actor.name}"` });
        continue;
      }
    }
    if (act.args.object !== undefined) {
      if (!table.object) {
        errors.push({ line: act.line, message: `"${act.verb}" doesn't take "${act.args.object.ref}" here; use key=value` });
      } else if (act.args[table.object] !== undefined) {
        errors.push({ line: act.line, message: `"${table.object}" is given twice` });
      } else {
        act.args[table.object] = act.args.object;
      }
      delete act.args.object;
    }
    checkBag(act.args, { ...table.args, ...EASE_ARG }, `"${act.verb}"`, act.line, errors);
    for (const [key, spec] of Object.entries(table.args)) {
      if (spec.required && act.args[key] === undefined && errors.length === errCount) {
        errors.push({ line: act.line, message: `"${act.verb}" needs ${key}=` });
      }
    }
    for (const r of collectRefs(act.args)) checkRef(r, act.line);
    if (act.t1 === null) act.t1 = act.t0 + table.defaultDuration;
  }

  errors.sort((a, b) => a.line - b.line);
  return { errors, warnings };
}
