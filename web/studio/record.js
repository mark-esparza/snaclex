// SnaCleX Studio — deterministic video export.
//
// Frames are rendered at a fixed time step (frame i shows t = i / fps), encoded
// with WebCodecs, and written into a WebM or MP4 container by a vendored muxer.
// Because each frame is rendered on demand rather than captured in real time,
// a slow machine produces the same smooth video as a fast one — it just takes
// longer.

import * as WebM from '../vendor/webm-muxer.js';
import * as MP4 from '../vendor/mp4-muxer.js';
import { stateAt } from './timeline.js';
import { drawOverlay } from './render/overlay.js';

export const EXPORT_SIZES = { '720p': [1280, 720], '1080p': [1920, 1080] };
export const FPS = 30;

// Codec candidates per container: [WebCodecs codec string, muxer codec id].
const FORMATS = {
  webm: { mime: 'video/webm', ext: 'webm', codecs: [['vp09.00.40.08', 'V_VP9'], ['vp8', 'V_VP8']] },
  mp4: { mime: 'video/mp4', ext: 'mp4', codecs: [['avc1.640028', 'avc']] },
};

export const formatInfo = (format) => FORMATS[format];

const bitrateFor = (width) => (width >= 1920 ? 12e6 : 6e6);

async function pickCodec(format, width, height, fps, VideoEncoder) {
  if (!VideoEncoder || typeof VideoEncoder.isConfigSupported !== 'function') return null;
  for (const [codec, muxCodec] of FORMATS[format].codecs) {
    try {
      const r = await VideoEncoder.isConfigSupported({ codec, width, height, bitrate: bitrateFor(width), framerate: fps });
      if (r && r.supported) return { codec, muxCodec };
    } catch { /* unsupported config shape: try the next codec */ }
  }
  return null;
}

// Which containers this browser can produce (checked at 720p).
export async function supportedFormats(env = globalThis) {
  const out = [];
  for (const format of Object.keys(FORMATS)) {
    if (await pickCodec(format, 1280, 720, FPS, env.VideoEncoder)) out.push(format);
  }
  return out;
}

const defaultDeps = {
  get VideoEncoder() { return globalThis.VideoEncoder; },
  get VideoFrame() { return globalThis.VideoFrame; },
  muxers: { webm: WebM, mp4: MP4 },
  stateAt,
  drawOverlay,
  createCanvas(w, h) {
    const c = document.createElement('canvas');
    c.width = w;
    c.height = h;
    return c;
  },
  yieldFn: () => new Promise((r) => setTimeout(r, 0)),
};

function abortError() {
  return new DOMException('Export cancelled', 'AbortError');
}

export async function exportVideo(opts, deps = defaultDeps) {
  const { timeline: tl, stage, width, height, format, fps = FPS, onProgress, signal, restore } = opts;
  let encoder = null;
  try {
    const picked = await pickCodec(format, width, height, fps, deps.VideoEncoder);
    if (!picked) throw new Error(`this browser cannot encode ${format.toUpperCase()} video at ${width}x${height}`);
    const { Muxer, ArrayBufferTarget } = deps.muxers[format];
    const target = new ArrayBufferTarget();
    const muxer = new Muxer({
      target,
      video: format === 'webm'
        ? { codec: picked.muxCodec, width, height, frameRate: fps }
        : { codec: picked.muxCodec, width, height },
      ...(format === 'mp4' ? { fastStart: 'in-memory' } : {}),
    });
    let encodeError = null;
    encoder = new deps.VideoEncoder({
      output: (chunk, meta) => muxer.addVideoChunk(chunk, meta),
      error: (e) => { encodeError = e; },
    });
    encoder.configure({ codec: picked.codec, width, height, bitrate: bitrateFor(width), framerate: fps });

    stage.setSize(width, height, 1);
    const overlay = deps.createCanvas(width, height);
    const octx = overlay.getContext('2d');
    const composite = deps.createCanvas(width, height);
    const cctx = composite.getContext('2d');

    const total = Math.ceil(tl.duration * fps - 1e-9);
    const frameUs = 1e6 / fps;
    for (let i = 0; i < total; i++) {
      if (signal && signal.aborted) throw abortError();
      if (encodeError) throw encodeError;
      const f = deps.stateAt(tl, i / fps);
      stage.apply(f);
      if (signal && signal.aborted) throw abortError();
      stage.render();
      deps.drawOverlay(octx, f, tl, stage, { width, height, style: tl.style });
      cctx.drawImage(stage.canvas, 0, 0, width, height);
      cctx.drawImage(overlay, 0, 0);
      const frame = new deps.VideoFrame(composite, { timestamp: Math.round(i * frameUs), duration: Math.round(frameUs) });
      encoder.encode(frame, { keyFrame: i % (2 * fps) === 0 });
      frame.close();
      if (i % 60 === 59) await encoder.flush();
      else if (i % 5 === 4) await deps.yieldFn();
      if (onProgress) onProgress((i + 1) / total);
    }
    await encoder.flush();
    if (encodeError) throw encodeError;
    muxer.finalize();
    return new Blob([target.buffer], { type: FORMATS[format].mime });
  } finally {
    if (encoder && encoder.state !== 'closed') {
      try { encoder.close(); } catch { /* already closed */ }
    }
    if (restore) restore();
  }
}
