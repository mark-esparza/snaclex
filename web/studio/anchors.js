// SnaCleX Studio — geometry conventions shared by the timeline and renderer.
//
// World units: 1 unit = 1 nm for real structures. A helix lies along its local
// +x axis, centred on its position; a chromosome's long axis is local y. Both
// the timeline (to place things "at dna@0.45") and the renderer (to draw them)
// use the functions here, so an anchor always lands where the geometry is.
//
// Pure module: no DOM, no WebGL.

export const RISE_PER_TURN = 3.4;   // nm per helical turn (10 bp × 0.34 nm)
export const KINETOCHORE_OFFSET = 1.2;

export const add = (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
export const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
export const mul = (a, k) => [a[0] * k, a[1] * k, a[2] * k];
export const len = (a) => Math.hypot(a[0], a[1], a[2]);
export const lerp = (a, b, x) => a + (b - a) * x;
export const lerpV = (a, b, x) => [lerp(a[0], b[0], x), lerp(a[1], b[1], x), lerp(a[2], b[2], x)];
export const normalize = (a) => { const l = len(a); return l > 1e-12 ? mul(a, 1 / l) : [0, 0, 1]; };

const RAD = Math.PI / 180;

// Rotate v by Euler angles in degrees, XYZ order (as three.js 'XYZ'):
// v' = Rx · Ry · Rz · v.
export function rotate(v, rotDeg) {
  let [x, y, z] = v;
  const [ax, ay, az] = rotDeg.map((d) => d * RAD);
  if (az) { const c = Math.cos(az), s = Math.sin(az); [x, y] = [c * x - s * y, s * x + c * y]; }
  if (ay) { const c = Math.cos(ay), s = Math.sin(ay); [x, z] = [c * x + s * z, -s * x + c * z]; }
  if (ax) { const c = Math.cos(ax), s = Math.sin(ax); [y, z] = [c * y - s * z, s * y + c * z]; }
  const clean = (n) => (Math.abs(n) < 1e-12 ? 0 : n);
  return [clean(x), clean(y), clean(z)];
}

export const helixLength = (turns) => turns * RISE_PER_TURN;

// Point and unit tangent at fraction u along a helix axis, in the helix's own
// frame (before position/rotation/scale). `curve` 0..1 bends the axis into an
// arc of up to 180°, keeping u = 0.5 at the origin with tangent +x.
export function helixFrame(u, turns, curve = 0) {
  const L = helixLength(turns);
  if (!curve) return { point: [(u - 0.5) * L, 0, 0], tangent: [1, 0, 0] };
  const angle = curve * Math.PI;
  const R = L / angle;
  const th = (u - 0.5) * angle;
  return { point: [R * Math.sin(th), R * (1 - Math.cos(th)), 0], tangent: [Math.cos(th), Math.sin(th), 0] };
}

// Local offset of a helix point, accounting for a cut that pushes the two
// halves apart along the axis.
function helixLocal(u, params, p) {
  const f = helixFrame(u, params.turns ?? 20, params.curve ?? 0);
  const cut = p.cut ?? 0;
  if (!cut) return f.point;
  const shift = ((p.gap ?? 4) * cut) / 2;
  return add(f.point, mul(f.tangent, u < (p.cutAt ?? 0.5) ? -shift : shift));
}

// Bounding radius of an actor, in world units.
export function actorExtent(type, params, stateScale) {
  const k = stateScale ?? params.scale ?? 1;
  switch (type) {
    case 'helix': return (helixLength(params.turns ?? 20) / 2) * k;
    case 'chromosome': return ((params.length ?? 12) / 2) * k;
    case 'cell': return (params.radius ?? 40) * k;
    case 'ring': return (params.radius ?? 3) * k;
    case 'protein': return ((params.size ?? 6) / 2) * k;
    case 'strand': return ((params.length ?? 10) / 2) * k;
    default: return 1 * k;
  }
}

// World position of a reference on an actor in a given state.
//   {ref}            -> the actor's position
//   {ref, u}         -> fraction u along a helix
//   {ref, part}      -> a named part (helix start/end, chromosome
//                       centromere/kinetochore/top/bottom, strand/fibre tip)
export function anchorPoint(decl, state, ref) {
  const p = state.p || {};
  let local = null;
  if (decl.type === 'helix' && (ref.u !== undefined || ref.part === 'start' || ref.part === 'end')) {
    const u = ref.u !== undefined ? ref.u : ref.part === 'start' ? 0 : 1;
    local = helixLocal(u, decl.params, p);
  } else if (ref.part === 'tip' && state.ends && state.ends.tip) {
    return [...state.ends.tip];
  } else if (decl.type === 'chromosome' && ref.part) {
    const L = (decl.params.length ?? 12) * (1.4 - 0.4 * (p.condense ?? 1));
    local = {
      centromere: [0, 0, 0],
      kinetochore: [-KINETOCHORE_OFFSET, 0, 0],
      top: [0, L / 2, 0],
      bottom: [0, -L / 2, 0],
    }[ref.part] || [0, 0, 0];
  }
  if (!local) return [...state.pos];
  return add(state.pos, rotate(mul(local, state.scale ?? 1), state.rot || [0, 0, 0]));
}
