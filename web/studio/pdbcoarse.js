// SnaCleX Studio — coarse-grain a protein chain into soft "beads".
//
// The animation look wants smooth blobs, not atoms. We take one Cα per
// residue (PDB fixed columns, first model, first alternate location), average
// each run of consecutive residues into a bead, and centre the result on the
// origin. Coordinates are converted from Å to nm (world units).
//
// Pure module: no DOM, no WebGL.

export const BEAD_RADIUS = 0.9;

export function caCoords(pdbText, chain) {
  const out = [];
  const seen = new Set();
  for (const line of String(pdbText).split(/\r?\n/)) {
    if (line.startsWith('ENDMDL')) break;
    if (!line.startsWith('ATOM  ')) continue;
    if (line.slice(12, 16).trim() !== 'CA') continue;
    if (line[21] !== chain) continue;
    const alt = line[16];
    if (alt !== ' ' && alt !== 'A' && alt !== undefined) continue;
    const resKey = line.slice(22, 27);
    if (seen.has(resKey)) continue;
    seen.add(resKey);
    const x = parseFloat(line.slice(30, 38));
    const y = parseFloat(line.slice(38, 46));
    const z = parseFloat(line.slice(46, 54));
    if ([x, y, z].every(Number.isFinite)) out.push([x / 10, y / 10, z / 10]);
  }
  return out;
}

export function coarseBeads(coords, perBead = 4) {
  const centres = [];
  for (let i = 0; i < coords.length; i += perBead) {
    const group = coords.slice(i, i + perBead);
    const c = [0, 0, 0];
    for (const p of group) { c[0] += p[0]; c[1] += p[1]; c[2] += p[2]; }
    centres.push(c.map((v) => v / group.length));
  }
  if (!centres.length) return { beads: [], radius: 0 };
  const mid = [0, 1, 2].map((k) => centres.reduce((s, c) => s + c[k], 0) / centres.length);
  let radius = 0;
  const beads = centres.map((c) => {
    const p = [c[0] - mid[0], c[1] - mid[1], c[2] - mid[2]];
    radius = Math.max(radius, Math.hypot(...p));
    return [p[0], p[1], p[2], BEAD_RADIUS];
  });
  return { beads, radius: radius + BEAD_RADIUS };
}
