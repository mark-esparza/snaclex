// SnaCleX Studio — the 2D text layer over the 3D preview: captions, actor
// labels with leader lines, the scene note ("Stylised; not to scale") and the
// closing credits card. Drawn on a separate canvas for the preview and
// composited onto each frame for export, so text is always crisp.
//
// layoutCaption() is pure (unit-tested); drawOverlay() needs a 2D context and
// a Stage for projecting world points to the screen.

import { actorExtent } from '../anchors.js';

const FONT = '"Segoe UI", "Helvetica Neue", Helvetica, Arial, sans-serif';

// Word-wrap text to lines no wider than width × maxWidthFrac. Newlines in the
// text force breaks; a single over-long word gets a line of its own.
export function layoutCaption(text, width, measure, maxWidthFrac = 0.8) {
  const max = width * maxWidthFrac;
  const lines = [];
  for (const para of String(text).split('\n')) {
    let line = '';
    for (const word of para.split(/\s+/).filter(Boolean)) {
      const next = line ? `${line} ${word}` : word;
      if (line && measure(next) > max) { lines.push(line); line = word; } else line = next;
    }
    lines.push(line);
  }
  return lines;
}

function roundRect(ctx, x, y, w, h, r) {
  ctx.beginPath();
  ctx.moveTo(x + r, y);
  ctx.arcTo(x + w, y, x + w, y + h, r);
  ctx.arcTo(x + w, y + h, x, y + h, r);
  ctx.arcTo(x, y + h, x, y, r);
  ctx.arcTo(x, y, x + w, y, r);
  ctx.closePath();
}

function palette(style) {
  return style === 'light'
    ? { text: '#1b1b1f', band: 'rgba(255,255,255,0.78)', line: 'rgba(30,30,40,0.6)', card: 'rgba(255,255,255,0.88)' }
    : { text: '#f4f4f8', band: 'rgba(0,0,0,0.55)', line: 'rgba(230,235,255,0.65)', card: 'rgba(6,8,16,0.82)' };
}

// How far above an actor's centre its label sits: the radius across the
// actor, not along it (a long helix is thin).
function labelRadius(stage, name, decl, s) {
  if (decl.type === 'helix') return 2.4 * s.scale;
  if (decl.type === 'strand') return 1 * s.scale;
  return stage.viewRadius(name) ?? actorExtent(decl.type, decl.params, s.scale);
}

function drawLabels(ctx, f, tl, stage, W, H, pal) {
  const size = Math.round(H * 0.024);
  ctx.font = `600 ${size}px ${FONT}`;
  ctx.textAlign = 'center';
  ctx.textBaseline = 'bottom';
  for (const name of tl.order) {
    const decl = tl.actors[name];
    const s = f.actors[name];
    const text = decl.type === 'label' ? decl.params.text : decl.params.label;
    if (!text || !s || !s.visible || s.opacity < 0.05) continue;
    const base = decl.type === 'fibre' || decl.type === 'arrow' ? s.ends.tip || s.ends.to : s.pos;
    ctx.globalAlpha = Math.min(1, s.opacity);
    if (decl.type === 'label') {
      const at = stage.project(base);
      if (!at) continue;
      ctx.fillStyle = pal.text;
      ctx.fillText(text, at[0], at[1]);
      continue;
    }
    const r = labelRadius(stage, name, decl, s);
    const edge = stage.project([base[0], base[1] + r, base[2]]);
    const at = stage.project([base[0], base[1] + r * 1.5 + 1, base[2]]);
    if (!edge || !at) continue;
    // Keep the label inside the frame even when its actor nearly fills it.
    at[1] = Math.max(at[1], size * 1.8);
    at[0] = Math.min(Math.max(at[0], size * 4), W - size * 4);
    ctx.strokeStyle = pal.line;
    ctx.lineWidth = Math.max(1, H / 720);
    ctx.beginPath();
    ctx.moveTo(edge[0], edge[1]);
    ctx.lineTo(at[0], at[1] + 2);
    ctx.stroke();
    const w = ctx.measureText(text).width + size * 0.8;
    ctx.fillStyle = pal.band;
    roundRect(ctx, at[0] - w / 2, at[1] - size * 1.35, w, size * 1.35, size * 0.3);
    ctx.fill();
    ctx.fillStyle = pal.text;
    ctx.fillText(text, at[0], at[1] - size * 0.18);
  }
  ctx.globalAlpha = 1;
}

function drawCaptions(ctx, f, W, H, pal) {
  if (!f.captions.length) return;
  const size = Math.round(H * 0.038);
  ctx.font = `500 ${size}px ${FONT}`;
  const lines = f.captions.flatMap((c) => layoutCaption(c, W, (s) => ctx.measureText(s).width));
  const lh = size * 1.3;
  const pad = size * 0.6;
  const bandH = lines.length * lh + pad * 2;
  const y0 = H - bandH - H * 0.04;
  const bandW = Math.max(...lines.map((l) => ctx.measureText(l).width)) + pad * 3;
  ctx.fillStyle = pal.band;
  roundRect(ctx, (W - bandW) / 2, y0, bandW, bandH, size * 0.35);
  ctx.fill();
  ctx.fillStyle = pal.text;
  ctx.textAlign = 'center';
  ctx.textBaseline = 'middle';
  lines.forEach((l, i) => ctx.fillText(l, W / 2, y0 + pad + lh * (i + 0.5)));
}

function drawNote(ctx, tl, H, pal) {
  if (!tl.note) return;
  const size = Math.round(H * 0.02);
  ctx.font = `400 ${size}px ${FONT}`;
  ctx.globalAlpha = 0.6;
  ctx.fillStyle = pal.text;
  ctx.textAlign = 'left';
  ctx.textBaseline = 'top';
  ctx.fillText(tl.note, size, size);
  ctx.globalAlpha = 1;
}

function drawCredits(ctx, tl, W, H, pal) {
  const title = tl.title || 'SnaCleX Studio';
  const tSize = Math.round(H * 0.045), bSize = Math.round(H * 0.024);
  const measure = (font) => (s) => { ctx.font = font; return ctx.measureText(s).width; };
  const tFont = `700 ${tSize}px ${FONT}`, bFont = `400 ${bSize}px ${FONT}`;
  const body = layoutCaption(tl.creditsText.replace(' · ', '\n'), W, measure(bFont), 0.6);
  const cardW = Math.min(W * 0.75, Math.max(measure(tFont)(title), ...body.map(measure(bFont))) + tSize * 2);
  const cardH = tSize * 2.2 + body.length * bSize * 1.45 + bSize;
  const x = (W - cardW) / 2, y = (H - cardH) / 2;
  ctx.fillStyle = pal.card;
  roundRect(ctx, x, y, cardW, cardH, tSize * 0.3);
  ctx.fill();
  ctx.fillStyle = pal.text;
  ctx.textAlign = 'center';
  ctx.textBaseline = 'top';
  ctx.font = tFont;
  ctx.fillText(title, W / 2, y + tSize * 0.7);
  ctx.font = bFont;
  body.forEach((l, i) => ctx.fillText(l, W / 2, y + tSize * 2.1 + i * bSize * 1.45));
}

export function drawOverlay(ctx, f, tl, stage, { width, height, style }) {
  const pal = palette(style);
  ctx.clearRect(0, 0, width, height);
  drawLabels(ctx, f, tl, stage, width, height, pal);
  drawNote(ctx, tl, height, pal);
  if (f.creditsVisible && tl.creditsText) drawCredits(ctx, tl, width, height, pal);
  else drawCaptions(ctx, f, width, height, pal);
}
