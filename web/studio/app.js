// SnaCleX Animation Studio — page wiring.
//
// Temporary development hook (replaced by the full page in a later task):
// /studio.html?dev=1 renders a fixed scene and exposes window.__dev.seek(t).
import { parse } from './parse.js';
import { compile, stateAt } from './timeline.js';
import { Stage } from './render/stage.js';
import { loadStructures } from './structures.js';

const DEV_SCENE = `
actor dna = helix turns=20
actor p = protein at=0,25,0
actor sep = protein from=7NJ1 chain=A at=-25,-20,0
actor sec = protein from=7NJ1 chain=B at=0,-20,0
0-3s p move to=dna@0.5
3-5s dna unwind from=0.45 to=0.55
5-7s camera zoom p distance=30
8s camera cut sep distance=40
`;

if (new URLSearchParams(location.search).has('dev')) {
  const r = parse(DEV_SCENE);
  const c = compile(r.scene);
  const canvas = document.getElementById('gl');
  const stage = new Stage(canvas, { style: c.timeline.style });
  const rect = canvas.getBoundingClientRect();
  stage.setSize(Math.round(rect.width), Math.round(rect.height), Math.min(2, devicePixelRatio));
  const { beads } = await loadStructures(r.scene);
  stage.build(c.timeline, beads);
  const seek = (t) => { stage.apply(stateAt(c.timeline, t)); stage.render(); };
  seek(0);
  window.__dev = { seek, timeline: c.timeline, stage };
}
