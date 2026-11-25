// SnaCleX Studio — timeline compiler and frame evaluator.
//
// compile(scene) turns a validated Scene into a Timeline: every action is
// expanded into tracks, each animating one property of one actor (or of the
// camera) from a value to a value over [t0, t1]. stateAt(timeline, t) is a
// pure function of t that returns a fresh FrameState, so scrubbing, replaying
// and recording always show the same frame for the same time.
//
// References inside actions ("move to=dna@0.45") resolve when the script is
// compiled, against where things are at that moment. Attachments ("attach
// to=q"), rings around things, labels on things and fibre endpoints resolve
// every frame instead, so they follow what they are attached to.
//
// Pure module: no DOM, no WebGL.

import { ACTOR_TYPES } from './schema.js';
import {
  add, sub, mul, len, lerp, lerpV, normalize, rotate,
  actorExtent, anchorPoint, helixFrame,
} from './anchors.js';

export const FOV = 35;
export const OFFSCREEN = 60;
export const CREDITS_SECONDS = 4;
export const PULL_FRACTION = 0.7;
export const INITIAL_GLOW = 0.4;
const WIDE_DIR = normalize([0, 0.25, 1]);
const DIR_VEC = {
  left: [-1, 0, 0], right: [1, 0, 0], top: [0, 1, 0],
  bottom: [0, -1, 0], front: [0, 0, 1], back: [0, 0, -1],
};
const AXIS_VEC = { x: [1, 0, 0], y: [0, 1, 0], z: [0, 0, 1] };
const BASE_PROPS = ['visible', 'opacity', 'pos', 'rot', 'scale', 'color', 'glow'];
const EPS = 1e-9;

// ---------------------------------------------------------------------------
// Easing and interpolation

export function ease(name, x) {
  const t = Math.min(1, Math.max(0, x));
  switch (name) {
    case 'linear': return t;
    case 'in': return t * t * t;
    case 'out': return 1 - (1 - t) ** 3;
    default: return t < 0.5 ? 4 * t * t * t : 1 - (-2 * t + 2) ** 3 / 2;
  }
}

const hexToRgb = (h) => [1, 3, 5].map((i) => parseInt(h.slice(i, i + 2), 16));
const rgbToHex = (c) => '#' + c.map((v) => Math.round(Math.min(255, Math.max(0, v))).toString(16).padStart(2, '0')).join('');

function interp(a, b, x) {
  if (typeof a === 'number') return lerp(a, b, x);
  if (Array.isArray(a)) return lerpV(a, b, x);
  if (typeof a === 'string') return rgbToHex(lerpV(hexToRgb(a), hexToRgb(b), x));
  return x < 1 ? a : b; // booleans step at the end
}

const clone = (v) => (Array.isArray(v) ? [...v] : v);
const fmt = (n) => String(Math.round(n * 1000) / 1000);

// ---------------------------------------------------------------------------
// Evaluation

function newCtx(tl, t) {
  return { tl, t, cache: Object.create(null), visiting: new Set() };
}

function initialProp(tl, name, prop) {
  if (name === 'camera') return tl.camInitial ? tl.camInitial[prop] : prop === 'pos' ? [0, 12.5, 50] : [0, 0, 0];
  const init = tl.initial[name];
  return prop.startsWith('p.') ? init.p[prop.slice(2)] : init[prop];
}

function trackValue(track, t, ctx) {
  if (track.dynamic) {
    const d = track.dynamic;
    const at = d.until !== undefined ? Math.min(t, d.until) : t;
    const c = at === ctx.t ? ctx : newCtx(ctx.tl, at);
    return add(refPoint(d.ref, c), d.offset);
  }
  if (t >= track.t1 || track.t1 - track.t0 < EPS) {
    return track.orbit ? orbitPoint(track.orbit, 1) : track.to;
  }
  const x = ease(track.ease, (t - track.t0) / (track.t1 - track.t0));
  if (track.orbit) return orbitPoint(track.orbit, x);
  return interp(track.from, track.to, x);
}

function orbitPoint(o, x) {
  return add(o.center, rotate(o.offset, [0, o.degrees * x, 0]));
}

function propAt(tl, name, prop, t, ctx) {
  const list = tl.index[`${name}|${prop}`];
  if (list) {
    for (let i = list.length - 1; i >= 0; i--) {
      const track = tl.tracks[list[i]];
      if (track.t0 <= t + EPS) return trackValue(track, t, ctx);
    }
  }
  return initialProp(tl, name, prop);
}

// World point of a Ref at the context's time.
function refPoint(ref, ctx) {
  const decl = ctx.tl.actors[ref.ref];
  return anchorPoint(decl, actorState(ref.ref, ctx), ref);
}

function resolvePoint(v, ctx) {
  return Array.isArray(v) ? [...v] : refPoint(v, ctx);
}

function axisOf(ref, ctx) {
  const decl = ctx.tl.actors[ref.ref];
  const s = actorState(ref.ref, ctx);
  if (decl.type === 'helix') {
    const u = ref.u ?? 0.5;
    return normalize(rotate(helixFrame(u, decl.params.turns, decl.params.curve).tangent, s.rot));
  }
  if (decl.type === 'chromosome') return normalize(rotate([0, 1, 0], s.rot));
  return normalize(rotate([1, 0, 0], s.rot));
}

function actorState(name, ctx) {
  if (ctx.cache[name]) return ctx.cache[name];
  const { tl, t } = ctx;
  const init = tl.initial[name];
  if (ctx.visiting.has(name)) return init; // defensive: compile rejects cycles
  ctx.visiting.add(name);
  const s = { p: {} };
  for (const prop of BASE_PROPS) s[prop] = clone(propAt(tl, name, prop, t, ctx));
  for (const k of Object.keys(init.p)) s.p[k] = propAt(tl, name, `p.${k}`, t, ctx);

  const decl = tl.actors[name];
  const P = decl.params;
  if (decl.type === 'fibre' || decl.type === 'arrow') {
    const from = resolvePoint(P.from, ctx);
    const to = resolvePoint(P.to, ctx);
    s.ends = { from, to, tip: lerpV(from, to, s.p.grow ?? 1) };
    s.pos = [...from];
  } else if (decl.type === 'strand') {
    if (P.along) {
      const h = actorState(P.along.ref, ctx);
      const u0 = h.p.unwind0, u1 = h.p.unwind1;
      const uTip = lerp(u0, u1, s.p.grow);
      const tip = anchorPoint(tl.actors[P.along.ref], h, { ref: P.along.ref, u: uTip });
      s.ends = { along: P.along.ref, u0, u1, uTip, tip };
    } else {
      s.ends = { tip: add(s.pos, rotate([(P.length ?? 10) * s.p.grow * s.scale, 0, 0], s.rot)) };
    }
  } else if (decl.type === 'ring' && P.around) {
    s.ends = { axis: axisOf(P.around, ctx) };
  }
  ctx.visiting.delete(name);
  ctx.cache[name] = s;
  return s;
}

export function stateAt(tl, t) {
  const ctx = newCtx(tl, t);
  const actors = {};
  for (const name of tl.order) {
    const s = actorState(name, ctx);
    actors[name] = {
      ...s,
      pos: [...s.pos], rot: [...s.rot], p: { ...s.p },
      ...(s.ends ? { ends: JSON.parse(JSON.stringify(s.ends)) } : {}),
    };
  }
  const camera = {
    pos: [...propAt(tl, 'camera', 'pos', t, ctx)],
    target: [...propAt(tl, 'camera', 'target', t, ctx)],
    fov: FOV,
  };
  const captions = tl.captions.filter((c) => t >= c.t0 - EPS && t < c.t1 - EPS).map((c) => c.text);
  const creditsVisible = !!tl.creditsText && t >= tl.duration - CREDITS_SECONDS - EPS;
  return { t, actors, camera, captions, creditsVisible };
}

// ---------------------------------------------------------------------------
// Compilation

function refNames(v) {
  if (v && typeof v === 'object' && !Array.isArray(v) && typeof v.ref === 'string') return [v.ref];
  return [];
}

// Every edge along which one actor's position depends on another's, for
// cycle detection: at=, around=, on=, along=, fibre/arrow ends, attach.
function dependencyEdges(scene) {
  const edges = new Map(scene.actors.map((a) => [a.name, []]));
  for (const a of scene.actors) {
    for (const key of ['at', 'around', 'on', 'along', 'from', 'to']) {
      for (const r of refNames(a.params[key])) edges.get(a.name).push({ to: r, line: a.line });
    }
  }
  for (const act of scene.actions) {
    if (act.verb === 'attach' && edges.has(act.target)) {
      edges.get(act.target).push({ to: act.args.to.ref, line: act.line });
    }
  }
  return edges;
}

function findCycle(edges) {
  const state = new Map();
  const stack = [];
  let found = null;
  const visit = (n) => {
    if (found) return;
    state.set(n, 1);
    stack.push(n);
    for (const e of edges.get(n) || []) {
      if (state.get(e.to) === 1) {
        const i = stack.indexOf(e.to);
        found = { path: [...stack.slice(i), e.to], line: e.line };
        return;
      }
      if (!state.has(e.to)) visit(e.to);
      if (found) return;
    }
    stack.pop();
    state.set(n, 2);
  };
  for (const n of edges.keys()) if (!state.has(n)) visit(n);
  return found;
}

function wideFraming(tl, state) {
  let names = tl.order.filter((n) => state.actors[n].visible && state.actors[n].opacity > 0.01);
  if (!names.length) names = tl.order;
  if (!names.length) return { target: [0, 0, 0], pos: mul(WIDE_DIR, 50) };
  const centroid = mul(names.reduce((acc, n) => add(acc, state.actors[n].pos), [0, 0, 0]), 1 / names.length);
  let R = 1;
  for (const n of names) {
    const a = tl.actors[n];
    R = Math.max(R, actorExtent(a.type, a.params, state.actors[n].scale) + len(sub(state.actors[n].pos, centroid)));
  }
  const d = (R / Math.tan((FOV / 2) * (Math.PI / 180))) * 1.15;
  return { target: centroid, pos: add(centroid, mul(WIDE_DIR, d)) };
}

export function compile(scene) {
  const errors = [];
  const warnings = [];

  const cycle = findCycle(dependencyEdges(scene));
  if (cycle) {
    errors.push({ line: cycle.line, message: `attach cycle: ${cycle.path.join(' → ')}` });
    return { timeline: null, errors, warnings };
  }

  const tl = {
    title: scene.title, size: [...scene.size], style: scene.style, note: scene.note,
    duration: 0, creditsText: null,
    actors: {}, order: [], initial: {}, tracks: [], index: {},
    camInitial: null, captions: [], warnings,
  };
  for (const a of scene.actors) {
    tl.actors[a.name] = { type: a.type, params: a.params, line: a.line };
    tl.order.push(a.name);
  }

  // Initial states. `at=` may reference another actor's initial position.
  const firstVerb = {};
  const hasVerb = {};
  for (const act of [...scene.actions].sort((x, y) => x.t0 - y.t0 || x.line - y.line)) {
    if (!tl.actors[act.target]) continue;
    if (!firstVerb[act.target]) firstVerb[act.target] = act.verb;
    (hasVerb[act.target] ||= new Set()).add(act.verb);
  }
  const resolving = new Set();
  const initialOf = (name) => {
    if (tl.initial[name]) return tl.initial[name];
    resolving.add(name);
    const { type, params } = tl.actors[name];
    const hiddenAtStart = params.hidden || ['enter', 'show'].includes(firstVerb[name]);
    const verbs = hasVerb[name] || new Set();
    const p = { ...ACTOR_TYPES[type].stateProps };
    if ('grow' in p && verbs.has('grow')) p.grow = 0;
    if ('condense' in p && verbs.has('condense')) p.condense = 0;
    let pos = [0, 0, 0];
    if (Array.isArray(params.at)) pos = [...params.at];
    else if (params.at && params.at.ref) {
      const other = initialOf(params.at.ref);
      pos = anchorPoint(tl.actors[params.at.ref], other, params.at);
    }
    const s = {
      visible: !hiddenAtStart, opacity: hiddenAtStart ? 0 : 1, pos,
      rot: [...params.rot], scale: params.scale, color: params.color, glow: INITIAL_GLOW, p,
    };
    resolving.delete(name);
    tl.initial[name] = s;
    return s;
  };
  for (const name of tl.order) initialOf(name);

  // Track bookkeeping. One overlap warning per action at most.
  const warned = new Set();
  const pushTrack = (track, act) => {
    const key = `${track.actor}|${track.prop}`;
    const list = (tl.index[key] ||= []);
    if (track.from === undefined && !track.dynamic && !track.orbit) {
      track.from = clone(propAt(tl, track.actor, track.prop, track.t0, newCtx(tl, track.t0)));
    }
    if (act && !warned.has(act) && track.t1 - track.t0 > EPS) {
      for (let i = list.length - 1; i >= 0; i--) {
        const prev = tl.tracks[list[i]];
        if (prev.line !== track.line && prev.t1 - prev.t0 > EPS && prev.t1 > track.t0 + EPS) {
          const prop = track.prop.startsWith('p.') ? track.prop.slice(2) : track.prop;
          const who = track.actor === 'camera' ? 'the camera' : `"${track.actor}"`;
          warnings.push({ line: track.line, message: `line ${track.line} overrides line ${prev.line} (${prop} of ${who})` });
          warned.add(act);
          break;
        }
      }
    }
    list.push(tl.tracks.length);
    tl.tracks.push(track);
  };
  const at = (t) => stateAt(tl, t);
  const point = (v, t) => (Array.isArray(v) ? [...v] : anchorPoint(tl.actors[v.ref], at(t).actors[v.ref], v));
  const cur = (name, prop, t) => clone(propAt(tl, name, prop, t, newCtx(tl, t)));

  // Implicit attachments: rings around something, labels on something.
  for (const name of tl.order) {
    const { type, params, line } = tl.actors[name];
    const ref = type === 'ring' ? params.around : type === 'label' ? params.on : null;
    if (ref) pushTrack({ actor: name, prop: 'pos', t0: 0, t1: 0, line, dynamic: { ref, offset: [0, 0, 0] } });
  }

  const sorted = [...scene.actions].sort((x, y) => x.t0 - y.t0 || x.line - y.line);
  for (const act of sorted.filter((a) => a.target !== 'camera')) {
    const { t0, t1, line, args } = act;
    const name = act.target;
    const e = args.ease || 'smooth';
    const tween = (prop, to, extra = {}) => pushTrack({ actor: name, prop, t0, t1, to, ease: e, line, ...extra }, act);
    const step = (prop, to, when = t0) => pushTrack({ actor: name, prop, t0: when, t1: when, from: to, to, ease: e, line });

    if (name === 'caption') { tl.captions.push({ t0, t1, text: act.text }); continue; }
    switch (act.verb) {
      case 'enter': {
        const to = args.to ? point(args.to, t1) : cur(name, 'pos', t0);
        const from = args.from.dir ? add(to, mul(DIR_VEC[args.from.dir], OFFSCREEN)) : point(args.from, t0);
        step('visible', true);
        tween('pos', to, { from });
        tween('opacity', 1, { from: 0 });
        break;
      }
      case 'leave': {
        const start = cur(name, 'pos', t0);
        const to = args.to.dir ? add(start, mul(DIR_VEC[args.to.dir], OFFSCREEN)) : point(args.to, t1);
        tween('pos', to);
        tween('opacity', 0);
        step('visible', false, t1);
        break;
      }
      case 'move': tween('pos', point(args.to, t1)); break;
      case 'attach': {
        const offset = sub(cur(name, 'pos', t0), point(args.to, t0));
        pushTrack({ actor: name, prop: 'pos', t0, t1: t0, line, dynamic: { ref: args.to, offset } }, act);
        break;
      }
      case 'detach': step('pos', cur(name, 'pos', t0)); break;
      case 'rotate': tween('rot', add(cur(name, 'rot', t0), mul(AXIS_VEC[args.axis], args.by))); break;
      case 'scale': tween('scale', args.to); break;
      case 'fade': tween('opacity', args.to); break;
      case 'glow': tween('glow', args.to); break;
      case 'color': tween('color', args.to); break;
      case 'show': step('visible', true); tween('opacity', 1); break;
      case 'hide': tween('opacity', 0); step('visible', false, t1); break;
      case 'unwind': {
        const c0 = cur(name, 'p.unwind0', t0), c1 = cur(name, 'p.unwind1', t0);
        const closed = Math.abs(c1 - c0) < EPS;
        const mid = (args.from + args.to) / 2;
        tween('p.unwind0', args.from, closed ? { from: mid } : {});
        tween('p.unwind1', args.to, closed ? { from: mid } : {});
        break;
      }
      case 'rewind': {
        const mid = (cur(name, 'p.unwind0', t0) + cur(name, 'p.unwind1', t0)) / 2;
        tween('p.unwind0', mid);
        tween('p.unwind1', mid);
        break;
      }
      case 'cut': step('p.cutAt', args.at); step('p.gap', args.gap); tween('p.cut', 1); break;
      case 'grow': tween('p.grow', args.to); break;
      case 'shrink': tween('p.grow', 0); break;
      case 'peel': tween('p.peel', 1); break;
      case 'open': tween('p.open', 1); break;
      case 'close': tween('p.open', 0); break;
      case 'cleave': tween('p.cleave', 1); break;
      case 'separate': step('p.gap', args.gap); tween('p.separate', 1); break;
      case 'condense': tween('p.condense', args.to); break;
      case 'pair': {
        const w = at(t1).actors[args.with.ref];
        tween('pos', add(w.pos, rotate([args.gap, 0, 0], w.rot)));
        break;
      }
      case 'pull': {
        const target = args.target.ref;
        const start = cur(target, 'pos', t0);
        const end = lerpV(start, point(args.to, t1), PULL_FRACTION);
        pushTrack({ actor: target, prop: 'pos', t0, t1, to: end, ease: e, line }, act);
        break;
      }
      case 'pinch': tween('p.pinch', args.to); break;
      case 'split': step('p.gap', args.gap); tween('p.split', 1); break;
      default: break;
    }
  }

  // Credits and duration.
  const credited = scene.actors.filter((a) => a.type === 'protein' && a.params.from);
  if (scene.credits && credited.length) {
    tl.creditsText = 'Structures: '
      + credited.map((a) => `${a.params.from} chain ${a.params.chain}${a.params.label ? ` (${a.params.label})` : ''}`).join(', ')
      + ' · Rendered with SnaCleX Studio';
  }
  const last = scene.actions.reduce((m, a) => (a.t1 > m.t1 ? a : m), { t1: 0, line: 0 });
  if (scene.duration !== null) {
    for (const a of scene.actions) {
      if (a.t1 > scene.duration + EPS) {
        errors.push({ line: a.line, message: `action ends at ${fmt(a.t1)}s but duration is ${fmt(scene.duration)}s` });
      }
    }
    tl.duration = scene.duration;
  } else {
    tl.duration = last.t1 + 1 + (tl.creditsText ? CREDITS_SECONDS : 0);
  }

  // Camera: starts framing everything visible at t = 0, then follows its own
  // actions. Compiled after the actors because it only ever looks at them.
  tl.camInitial = wideFraming(tl, at(0));
  for (const act of sorted.filter((a) => a.target === 'camera')) {
    const { t0, t1, line, args } = act;
    const e = args.ease || 'smooth';
    const cam = (prop, extra) => pushTrack({ actor: 'camera', prop, t0, t1, ease: e, line, ...extra }, act);
    const curPos = cur('camera', 'pos', t0);
    const curTarget = cur('camera', 'target', t0);
    const dir = normalize(sub(curPos, curTarget));
    const distFor = (ref, t, given) => {
      if (given !== undefined) return given;
      const a = tl.actors[ref.ref];
      return 3 * actorExtent(a.type, a.params, at(t).actors[ref.ref].scale);
    };
    switch (act.verb) {
      case 'wide': {
        const f = wideFraming(tl, at(t1));
        cam('target', { to: f.target });
        cam('pos', { to: f.pos });
        break;
      }
      case 'zoom': {
        const tp = point(args.to, t1);
        cam('target', { to: tp });
        cam('pos', { to: add(tp, mul(dir, distFor(args.to, t1, args.distance))) });
        break;
      }
      case 'orbit': {
        const center = args.around ? point(args.around, t0) : curTarget;
        const offset = sub(curPos, center);
        cam('target', { to: center });
        cam('pos', { orbit: { center, offset, degrees: args.degrees }, to: orbitPoint({ center, offset, degrees: args.degrees }, 1) });
        break;
      }
      case 'follow': {
        const anchor0 = point(args.ref, t0);
        cam('target', { dynamic: { ref: args.ref, offset: [0, 0, 0], until: t1 } });
        cam('pos', { dynamic: { ref: args.ref, offset: sub(curPos, anchor0), until: t1 } });
        break;
      }
      case 'cut': {
        const tp = point(args.to, t0);
        const d = distFor(args.to, t0, args.distance);
        pushTrack({ actor: 'camera', prop: 'target', t0, t1: t0, from: tp, to: tp, line });
        pushTrack({ actor: 'camera', prop: 'pos', t0, t1: t0, from: add(tp, mul(dir, d)), to: add(tp, mul(dir, d)), line });
        break;
      }
      default: break;
    }
  }

  warnings.sort((a, b) => a.line - b.line);
  if (errors.length) return { timeline: null, errors, warnings };
  return { timeline: tl, errors, warnings };
}
