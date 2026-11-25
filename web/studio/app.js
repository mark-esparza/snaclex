// SnaCleX Animation Studio — page wiring.
//
// Script text -> parse -> fetch real structures -> compile -> Stage. Editing
// re-runs that pipeline (debounced); if the new text has errors, the last good
// scene keeps playing and the errors are listed by line. Playback, scrubbing
// and export all go through stateAt(timeline, t), so every view of time t is
// the same frame.

import { parse } from './parse.js';
import { compile, stateAt } from './timeline.js';
import { loadStructures } from './structures.js';
import { Stage } from './render/stage.js';
import { drawOverlay } from './render/overlay.js';
import { initExport } from './export-ui.js';

const $ = (id) => document.getElementById(id);
const els = {
  editor: $('editor'), diagnostics: $('diagnostics'), preview: $('preview'),
  gl: $('gl'), overlay: $('overlay'), play: $('play'), scrub: $('scrub'),
  clock: $('clock'), examples: $('examples'), open: $('open'), save: $('save'),
};

const EXAMPLES_DIR = 'studio/examples/';
const DEBOUNCE_MS = 300;

const state = {
  tl: null,
  beads: new Map(),
  t: 0,
  playing: false,
  dirty: false,         // editor changed since last load/save
  version: 0,           // bumps on every edit; stale async rebuilds are dropped
  builtVersion: -1,     // version of the text the current diagnostics describe
  needsRender: true,
  css: { w: 0, h: 0 },
  busy: false,          // true while exporting: editing and playback paused
};

const stage = new Stage(els.gl, { style: 'dark' });
const octx = els.overlay.getContext('2d');

// ---------------------------------------------------------------------------
// Sizing and rendering

function resize() {
  if (state.busy) return;
  const r = els.preview.getBoundingClientRect();
  const w = Math.max(1, Math.round(r.width)), h = Math.max(1, Math.round(r.height));
  const dpr = Math.min(2, window.devicePixelRatio || 1);
  state.css = { w, h };
  stage.setSize(w, h, dpr);
  els.overlay.width = Math.round(w * dpr);
  els.overlay.height = Math.round(h * dpr);
  octx.setTransform(dpr, 0, 0, dpr, 0, 0);
  state.needsRender = true;
}

export function renderAt(t) {
  if (!state.tl) {
    octx.clearRect(0, 0, state.css.w, state.css.h);
    stage.render();
    return;
  }
  const f = stateAt(state.tl, t);
  stage.apply(f);
  stage.render();
  drawOverlay(octx, f, state.tl, stage, { width: state.css.w, height: state.css.h, style: state.tl.style });
}

const fmtClock = (t, decimals) => {
  const m = Math.floor(t / 60);
  const s = t - m * 60;
  return `${m}:${decimals ? s.toFixed(1).padStart(4, '0') : String(Math.floor(s)).padStart(2, '0')}`;
};

function updateClock() {
  const d = state.tl ? state.tl.duration : 0;
  els.clock.textContent = `${fmtClock(state.t, true)} / ${fmtClock(d, false)}`;
  els.scrub.max = String(d || 1);
  els.scrub.value = String(state.t);
  els.play.textContent = state.playing ? '❚❚' : '▶';
  els.play.setAttribute('aria-label', state.playing ? 'Pause' : 'Play');
}

let lastFrame = performance.now();
function tick(now) {
  const dt = Math.min(0.1, (now - lastFrame) / 1000);
  lastFrame = now;
  if (state.playing && state.tl && !state.busy) {
    state.t += dt;
    if (state.t >= state.tl.duration) { state.t = state.tl.duration; state.playing = false; }
    state.needsRender = true;
  }
  if (state.needsRender && !state.busy) {
    state.needsRender = false;
    renderAt(state.t);
    updateClock();
  }
  requestAnimationFrame(tick);
}

function seek(t) {
  state.t = Math.min(Math.max(0, t), state.tl ? state.tl.duration : 0);
  state.needsRender = true;
  updateClock();
}

// ---------------------------------------------------------------------------
// Diagnostics

function lineRange(text, line) {
  const lines = text.split('\n');
  let start = 0;
  for (let i = 0; i < Math.min(line - 1, lines.length - 1); i++) start += lines[i].length + 1;
  return [start, start + (lines[line - 1] || '').length];
}

function showDiagnostics(errors, warnings) {
  els.diagnostics.replaceChildren();
  const add = (d, cls) => {
    const el = document.createElement('span');
    el.className = cls;
    el.textContent = `Line ${d.line}: ${d.message}`;
    el.addEventListener('click', () => {
      const [a, b] = lineRange(els.editor.value, d.line);
      els.editor.focus();
      els.editor.setSelectionRange(a, b);
    });
    els.diagnostics.appendChild(el);
  };
  errors.forEach((d) => add(d, 'err'));
  warnings.forEach((d) => add(d, 'warn'));
}

// ---------------------------------------------------------------------------
// Rebuilding the scene from the script

async function rebuild() {
  const version = ++state.version;
  const text = els.editor.value;
  const p = parse(text);
  if (p.errors.length) {
    showDiagnostics(p.errors, p.warnings);
    state.builtVersion = version;
    return;
  }
  if (!p.scene.actors.length && !p.scene.actions.length) {
    showDiagnostics([], []);
    state.builtVersion = version;
    return;
  }
  const { beads, warnings: structWarnings } = await loadStructures(p.scene);
  if (version !== state.version) return; // a newer edit superseded this one
  // A real structure knows its own size; framing and labels should use it.
  for (const a of p.scene.actors) {
    const b = beads.get(a.name);
    if (a.type === 'protein' && b) a.params.size = 2 * b.radius;
  }
  const c = compile(p.scene);
  if (c.errors.length) {
    showDiagnostics(c.errors, [...p.warnings, ...structWarnings, ...c.warnings]);
    state.builtVersion = version;
    return;
  }
  state.tl = c.timeline;
  state.beads = beads;
  stage.build(state.tl, beads);
  els.preview.style.background = state.tl.style === 'light' ? '#f4f1ea' : '#05060a';
  state.t = Math.min(state.t, state.tl.duration);
  showDiagnostics([], [...p.warnings, ...structWarnings, ...c.warnings].sort((a, b) => a.line - b.line));
  state.builtVersion = version;
  state.needsRender = true;
  updateClock();
}

let debounce = null;
function scheduleRebuild(delay = DEBOUNCE_MS) {
  clearTimeout(debounce);
  state.version++; // anything in flight is now stale
  debounce = setTimeout(rebuild, delay);
}

// ---------------------------------------------------------------------------
// Toolbar: examples, open, save

function slug(s) {
  return String(s || '').toLowerCase().replace(/[^a-z0-9]+/g, '-').replace(/^-+|-+$/g, '') || 'scene';
}

function setScript(text) {
  els.editor.value = text;
  state.dirty = false;
  state.t = 0;
  state.playing = false;
  scheduleRebuild(0);
}

els.editor.addEventListener('input', () => { state.dirty = true; scheduleRebuild(); });

els.examples.addEventListener('change', async () => {
  const name = els.examples.value;
  if (!name) return;
  if (state.dirty && els.editor.value.trim() && !window.confirm('Replace your unsaved script with this example?')) {
    els.examples.value = '';
    return;
  }
  try {
    const res = await fetch(`${EXAMPLES_DIR}${name}.txt`);
    if (!res.ok) throw new Error(`HTTP ${res.status}`);
    setScript(await res.text());
  } catch (e) {
    showDiagnostics([{ line: 1, message: `could not load the example (${e.message})` }], []);
  }
  els.examples.value = '';
});

els.open.addEventListener('change', async () => {
  const file = els.open.files && els.open.files[0];
  if (!file) return;
  setScript(await file.text());
  els.open.value = '';
});

function download(blob, filename) {
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url;
  a.download = filename;
  document.body.appendChild(a);
  a.click();
  a.remove();
  setTimeout(() => URL.revokeObjectURL(url), 10000);
}

els.save.addEventListener('click', () => {
  const title = state.tl ? state.tl.title : '';
  download(new Blob([els.editor.value], { type: 'text/plain' }), `${slug(title)}.txt`);
  state.dirty = false;
});

// ---------------------------------------------------------------------------
// Playback controls

els.play.addEventListener('click', () => {
  if (!state.tl) return;
  if (!state.playing && state.t >= state.tl.duration - 1e-6) state.t = 0;
  state.playing = !state.playing;
  lastFrame = performance.now();
  updateClock();
});
els.scrub.addEventListener('input', () => {
  state.playing = false;
  seek(parseFloat(els.scrub.value));
});

// ---------------------------------------------------------------------------
// Start

new ResizeObserver(resize).observe(els.preview);
resize();
updateClock();
requestAnimationFrame(tick);

initExport({
  stage, state, slug, download, renderAt,
  setBusy(busy) {
    state.busy = busy;
    state.playing = false;
    // Nothing may replace the scene while it is being exported.
    for (const el of [els.editor, els.examples, els.open]) el.disabled = busy;
    if (!busy) { resize(); state.needsRender = true; }
  },
});

if (new URLSearchParams(location.search).has('test')) {
  window.__studio = {
    getTimeline: () => (state.tl ? { order: [...state.tl.order], duration: state.tl.duration } : null),
    seek: (t) => { state.playing = false; seek(t); renderAt(state.t); },
    isPlaying: () => state.playing,
    time: () => state.t,
    ready: () => state.builtVersion === state.version && !!state.tl,
    structureCount: () => [...state.beads.values()].filter(Boolean).length,
  };
}
