import test from 'node:test';
import assert from 'node:assert/strict';
import { supportedFormats, exportVideo, EXPORT_SIZES } from '../record.js';

test('export sizes are exactly 720p and 1080p', () => {
  assert.deepEqual(EXPORT_SIZES, { '720p': [1280, 720], '1080p': [1920, 1080] });
});

test('supportedFormats with no VideoEncoder', async () => {
  assert.deepEqual(await supportedFormats({}), []);
});

test('supportedFormats filters by codec', async () => {
  const env = { VideoEncoder: { isConfigSupported: async (c) => ({ supported: c.codec.startsWith('vp') }) } };
  assert.deepEqual(await supportedFormats(env), ['webm']);
  const both = { VideoEncoder: { isConfigSupported: async () => ({ supported: true }) } };
  assert.deepEqual(await supportedFormats(both), ['webm', 'mp4']);
  const throws = { VideoEncoder: { isConfigSupported: async () => { throw new TypeError('bad'); } } };
  assert.deepEqual(await supportedFormats(throws), []);
});

// A fake world: records what the exporter does without a browser.
function fakeDeps(log) {
  class FakeEncoder {
    static async isConfigSupported() { return { supported: true }; }
    constructor({ output, error }) { this.output = output; this.error = error; log.encoder = this; }
    configure(cfg) { log.config = cfg; }
    encode(frame, opts) {
      log.frames.push({ timestamp: frame.timestamp, keyFrame: !!(opts && opts.keyFrame) });
      this.output({ byteLength: 10 }, {});
    }
    async flush() { log.flushes++; }
    close() { log.closed = true; }
  }
  class FakeFrame {
    constructor(source, init) { this.timestamp = init.timestamp; this.duration = init.duration; }
    close() { log.framesClosed++; }
  }
  class Target { constructor() { this.buffer = new ArrayBuffer(8); } }
  class Muxer {
    constructor(opts) { log.muxer = opts; }
    addVideoChunk() { log.chunks++; }
    finalize() { log.finalized = true; }
  }
  const ctx = { drawImage() {}, clearRect() {}, setTransform() {}, fillRect() {} };
  return {
    VideoEncoder: FakeEncoder,
    VideoFrame: FakeFrame,
    muxers: { webm: { Muxer, ArrayBufferTarget: Target }, mp4: { Muxer, ArrayBufferTarget: Target } },
    stateAt: (tl, t) => ({ t, actors: {}, camera: {}, captions: [], creditsVisible: false }),
    drawOverlay: () => {},
    createCanvas: (w, h) => ({ width: w, height: h, getContext: () => ctx }),
    yieldFn: async () => {},
  };
}

function fakeStage(log) {
  return {
    canvas: {},
    setSize(w, h, r) { log.sizes.push([w, h, r]); },
    apply(f) { log.applied.push(f.t); },
    render() {},
  };
}

const newLog = () => ({ frames: [], applied: [], sizes: [], flushes: 0, framesClosed: 0, chunks: 0 });

test('frame timestamps are fixed-step', async () => {
  const log = newLog();
  const blob = await exportVideo({
    timeline: { duration: 1, style: 'dark' }, stage: fakeStage(log),
    width: 1280, height: 720, format: 'webm', fps: 30,
  }, fakeDeps(log));
  assert.equal(log.frames.length, 30);
  assert.deepEqual(log.frames.slice(0, 3).map((f) => f.timestamp), [0, 33333, 66667]);
  assert.deepEqual(log.frames.map((f) => f.keyFrame), log.frames.map((_, i) => i % 60 === 0));
  assert.deepEqual(log.applied.slice(0, 2), [0, 1 / 30]);
  assert.equal(log.framesClosed, 30);
  assert.equal(log.config.width, 1280);
  assert.equal(log.config.bitrate, 6e6);
  assert.equal(log.config.framerate, 30);
  assert.equal(log.muxer.video.codec, 'V_VP9');
  assert.equal(log.finalized, true);
  assert.equal(blob.type, 'video/webm');
  assert.deepEqual(log.sizes[0], [1280, 720, 1]);
});

test('1080p mp4 uses avc, a higher bitrate, and fast start', async () => {
  const log = newLog();
  const blob = await exportVideo({
    timeline: { duration: 0.1, style: 'dark' }, stage: fakeStage(log),
    width: 1920, height: 1080, format: 'mp4',
  }, fakeDeps(log));
  assert.equal(log.config.bitrate, 12e6);
  assert.equal(log.config.codec, 'avc1.640028');
  assert.equal(log.muxer.video.codec, 'avc');
  assert.equal(log.muxer.fastStart, 'in-memory');
  assert.equal(blob.type, 'video/mp4');
});

test('keyframe every two seconds over a long export', async () => {
  const log = newLog();
  await exportVideo({ timeline: { duration: 5, style: 'dark' }, stage: fakeStage(log), width: 1280, height: 720, format: 'webm', fps: 30 }, fakeDeps(log));
  const keys = log.frames.map((f, i) => (f.keyFrame ? i : -1)).filter((i) => i >= 0);
  assert.deepEqual(keys, [0, 60, 120]);
  assert.ok(log.flushes >= 3, 'flushes periodically to bound memory');
});

test('abort stops export', async () => {
  const log = newLog();
  const ac = new AbortController();
  const deps = fakeDeps(log);
  const stage = fakeStage(log);
  const apply = stage.apply;
  stage.apply = (f) => { apply(f); if (log.applied.length === 5) ac.abort(); };
  await assert.rejects(
    exportVideo({ timeline: { duration: 2, style: 'dark' }, stage, width: 1280, height: 720, format: 'webm', signal: ac.signal }, deps),
    (e) => e.name === 'AbortError',
  );
  assert.ok(log.frames.length <= 6);
  assert.equal(log.closed, true);
  assert.equal(log.finalized, undefined);
});

test('restore runs even when export fails', async () => {
  const log = newLog();
  let restored = false;
  const deps = fakeDeps(log);
  deps.VideoFrame = class { constructor() { throw new Error('boom'); } };
  await assert.rejects(exportVideo({
    timeline: { duration: 1, style: 'dark' }, stage: fakeStage(log), width: 1280, height: 720,
    format: 'webm', restore: () => { restored = true; },
  }, deps), /boom/);
  assert.equal(restored, true);
});
