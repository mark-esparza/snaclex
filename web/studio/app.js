// SnaCleX Animation Studio — page wiring.
//
// Temporary development hook (replaced by the full page in a later task):
// /studio.html?dev=1 renders a fixed scene and exposes window.__dev.seek(t).
import { parse } from './parse.js';
import { compile, stateAt } from './timeline.js';
import { Stage } from './render/stage.js';
import { loadStructures } from './structures.js';

const DEV_SCENES = {
  '1': `
actor dna = helix turns=20
actor p = protein at=0,25,0
actor sep = protein from=7NJ1 chain=A at=-25,-20,0
actor sec = protein from=7NJ1 chain=B at=0,-20,0
0-3s p move to=dna@0.5
3-5s dna unwind from=0.45 to=0.55
5-7s camera zoom p distance=30
8s camera cut sep distance=40
`,
  '2': `
actor cell = cell radius=40
actor c1 = chromosome color=violet at=0,0,0
actor c2 = chromosome color=pink at=6,0,0
actor coh = ring around=c1 radius=2.6
actor pole = point at=-30,0,0 show=yes
actor f = fibre from=pole to=c1.kinetochore
0-1s f grow
1-2s coh open
2-3s coh cleave
2-3s c2 separate
3-5s cell pinch
6-8s cell split
`,
  '3': `
actor c1 = chromosome color=violet
actor coh = ring around=c1.top radius=2.6
actor dna = helix turns=6 at=0,-14,0
actor coh2 = ring around=dna@0.5
0s camera zoom c1 distance=34
1-2s coh open
1-2s coh2 open
3-4s coh cleave
3-4s c1 separate
`,
};

if (new URLSearchParams(location.search).has('dev')) {
  const r = parse(DEV_SCENES[new URLSearchParams(location.search).get('dev')] || DEV_SCENES['1']);
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
