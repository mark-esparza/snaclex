// SnaCleX Studio — views for molecular-scale actors: DNA helix, RNA strand,
// protein blob, invisible/visible point, arrow and (overlay-only) label.
//
// Each maker returns an ActorView {object, update(state, frame, timeline),
// dispose()}. update() is called every frame with the actor's ActorState and
// must be cheap when nothing changed: geometry is rebuilt only when the
// shape-changing state (unwind range, cut, growth) actually moves.

import * as THREE from '../../vendor/three.module.min.js';
import { MarchingCubes } from '../../vendor/MarchingCubes.js';
import { glowMaterial, setLook } from './materials.js';
import { helixFrame, rotate, add, mul, sub, len } from '../anchors.js';

const RAD = Math.PI / 180;
const changed = (a, b) => a === undefined || Math.abs(a - b) > 1e-3;

export function applyTransform(object, s) {
  object.position.set(s.pos[0], s.pos[1], s.pos[2]);
  object.rotation.set(s.rot[0] * RAD, s.rot[1] * RAD, s.rot[2] * RAD, 'XYZ');
  object.scale.setScalar(s.scale);
  object.visible = s.visible && s.opacity > 0.002;
}

function disposeTree(object) {
  object.traverse((o) => {
    if (o.geometry) o.geometry.dispose();
  });
}

function lighten(hex, k) {
  const c = new THREE.Color(hex);
  c.lerp(new THREE.Color('#ffffff'), k);
  return '#' + c.getHexString();
}

// ---------------------------------------------------------------------------
// Helix

const STRAND_RADIUS = 0.18;
const HELIX_RADIUS = 1;
const BUBBLE_RADIUS = 2.4;
const GROOVE_PHASE = 0.78 * Math.PI; // offset between strands: major/minor groove

function smoothstep(e0, e1, x) {
  const t = Math.min(1, Math.max(0, (x - e0) / (e1 - e0)));
  return t * t * (3 - 2 * t);
}

// 0 outside the unwound range, 1 inside, smooth at the edges.
function bubble(u, u0, u1, edge) {
  if (u1 - u0 < 1e-4) return 0;
  return smoothstep(u0 - edge, u0 + edge * 0.5, u) * (1 - smoothstep(u1 - edge * 0.5, u1 + edge, u));
}

// Point on strand k (0|1) at fraction u, in the helix's local frame.
export function strandPoint(params, p, k, u) {
  const turns = params.turns ?? 20;
  const f = helixFrame(u, turns, params.curve ?? 0);
  const t = f.tangent;
  const n1 = [-t[1], t[0], 0];
  const n2 = [0, 0, 1];
  const edge = 1.2 / turns; // ~1 turn of soft edge
  const b = bubble(u, p.unwind0 ?? 0, p.unwind1 ?? 0, edge);
  const r = HELIX_RADIUS + (BUBBLE_RADIUS - HELIX_RADIUS) * b;
  const phi = 2 * Math.PI * turns * u + (k ? GROOVE_PHASE : 0);
  let pt = add(f.point, add(mul(n1, r * Math.cos(phi)), mul(n2, r * Math.sin(phi))));
  // In the bubble the strands also part sideways so the opening reads clearly.
  pt = add(pt, mul(n2, (k ? -1 : 1) * 0.9 * b));
  const cut = p.cut ?? 0;
  if (cut) pt = add(pt, mul(t, (u < (p.cutAt ?? 0.5) ? -1 : 1) * ((p.gap ?? 4) * cut) / 2));
  return { pt, b };
}

function tube(points, radius, material, radial = 8) {
  const curve = new THREE.CatmullRomCurve3(points.map((q) => new THREE.Vector3(q[0], q[1], q[2])));
  const geo = new THREE.TubeGeometry(curve, Math.max(8, points.length * 2), radius, radial, false);
  return new THREE.Mesh(geo, material);
}

export function makeHelix(decl, _beads, style) {
  const group = new THREE.Group();
  const backbone = glowMaterial(decl.params.color, style);
  const bases = glowMaterial(lighten(decl.params.color, 0.45), style);
  const geomGroup = new THREE.Group();
  group.add(geomGroup);
  let last = {};

  function rebuild(p) {
    disposeTree(geomGroup);
    geomGroup.clear();
    const turns = decl.params.turns ?? 20;
    const samples = Math.ceil(turns * 20);
    const cut = (p.cut ?? 0) > 1e-3;
    const cutAt = p.cutAt ?? 0.5;
    for (const k of [0, 1]) {
      let seg = [];
      for (let i = 0; i <= samples; i++) {
        const u = i / samples;
        if (cut && seg.length && u >= cutAt && (i - 1) / samples < cutAt) {
          geomGroup.add(tube(seg, STRAND_RADIUS, backbone));
          seg = [];
        }
        seg.push(strandPoint(decl.params, p, k, u).pt);
      }
      if (seg.length > 1) geomGroup.add(tube(seg, STRAND_RADIUS, backbone));
    }
    // Base pairs: one rung per 0.34 nm; inside the bubble each strand keeps a
    // short stub (an unpaired base) instead of a full rung.
    const count = Math.round(turns * 10);
    const rung = new THREE.CylinderGeometry(0.11, 0.11, 1, 6);
    const mesh = new THREE.InstancedMesh(rung, bases, count * 2);
    const m = new THREE.Matrix4();
    const up = new THREE.Vector3(0, 1, 0);
    let n = 0;
    const place = (a, b) => {
      const va = new THREE.Vector3(...a), vb = new THREE.Vector3(...b);
      const dir = vb.clone().sub(va);
      const length = dir.length();
      if (length < 1e-4) return;
      const q = new THREE.Quaternion().setFromUnitVectors(up, dir.normalize());
      m.compose(va.clone().add(vb).multiplyScalar(0.5), q, new THREE.Vector3(1, length, 1));
      mesh.setMatrixAt(n++, m);
    };
    for (let i = 0; i < count; i++) {
      const u = (i + 0.5) / count;
      if (cut && Math.abs(u - cutAt) < 0.5 / count) continue;
      const a = strandPoint(decl.params, p, 0, u);
      const b = strandPoint(decl.params, p, 1, u);
      if (a.b < 0.3) {
        place(a.pt, b.pt);
      } else {
        const mid = mul(add(a.pt, b.pt), 0.5);
        const stub = (from) => add(from, mul(sub(mid, from), Math.min(1, 0.55 / Math.max(len(sub(mid, from)), 1e-6))));
        place(a.pt, stub(a.pt));
        place(b.pt, stub(b.pt));
      }
    }
    mesh.count = n;
    mesh.instanceMatrix.needsUpdate = true;
    geomGroup.add(mesh);
  }

  return {
    object: group,
    update(s) {
      applyTransform(group, s);
      setLook(backbone, s.color, s.opacity, s.glow);
      setLook(bases, lighten(s.color, 0.45), s.opacity, s.glow);
      const p = s.p;
      if (['unwind0', 'unwind1', 'cut', 'cutAt', 'gap'].some((k) => changed(last[k], p[k] ?? 0))) {
        rebuild(p);
        last = { ...p };
      }
    },
    dispose() { disposeTree(group); backbone.dispose(); bases.dispose(); },
  };
}

// ---------------------------------------------------------------------------
// Strand (RNA or single-stranded DNA)

export function makeStrand(decl, _beads, style) {
  const group = new THREE.Group();
  const mat = glowMaterial(decl.params.color, style);
  let key = '';
  let mesh = null;
  const rebuild = (pts) => {
    if (mesh) { group.remove(mesh); mesh.geometry.dispose(); mesh = null; }
    if (pts.length > 1) { mesh = tube(pts, 0.16, mat, 6); group.add(mesh); }
  };
  return {
    object: group,
    update(s, f, tl) {
      setLook(mat, s.color, s.opacity, s.glow);
      group.visible = s.visible && s.opacity > 0.002;
      const e = s.ends || {};
      if (e.along) {
        // World-space: follow the helix bubble, riding just outside it.
        group.position.set(0, 0, 0);
        group.rotation.set(0, 0, 0);
        group.scale.setScalar(1);
        const h = f.actors[e.along];
        const hd = tl.actors[e.along];
        const k = `${e.u0.toFixed(4)}|${e.uTip.toFixed(4)}|${h.pos.join()}|${h.rot.join()}|${h.scale}|${s.p.peel}`;
        if (k === key) return;
        key = k;
        const n = Math.max(2, Math.ceil(Math.abs(e.uTip - e.u0) * (hd.params.turns ?? 20) * 12));
        const pts = [];
        for (let i = 0; i <= n; i++) {
          const u = e.u0 + ((e.uTip - e.u0) * i) / n;
          const local = strandPoint(hd.params, h.p, 0, u).pt;
          const axisLocal = helixFrame(u, hd.params.turns ?? 20, hd.params.curve ?? 0).point;
          const out = mul(sub(local, axisLocal), 1.35 + (s.p.peel ?? 0) * 2);
          const world = add(h.pos, rotate(mul(add(axisLocal, out), h.scale), h.rot));
          pts.push(world);
        }
        rebuild(pts);
      } else {
        applyTransform(group, s);
        const L = (decl.params.length ?? 10) * s.p.grow;
        const k = `${L.toFixed(3)}`;
        if (k === key) return;
        key = k;
        const n = Math.max(2, Math.ceil(L * 3));
        const pts = [];
        for (let i = 0; i <= n; i++) {
          const x = (L * i) / n;
          pts.push([x, 0.35 * Math.sin(x * 1.3), 0.35 * Math.cos(x * 0.9)]);
        }
        rebuild(pts);
      }
    },
    dispose() { if (mesh) mesh.geometry.dispose(); mat.dispose(); },
  };
}

// ---------------------------------------------------------------------------
// Protein

// FNV-1a 32-bit hash and a small seeded PRNG, so generic blobs are stable.
export function fnv1a(str) {
  let h = 0x811c9dc5;
  for (let i = 0; i < str.length; i++) { h ^= str.charCodeAt(i); h = Math.imul(h, 0x01000193); }
  return h >>> 0;
}
function mulberry32(a) {
  return () => {
    a |= 0; a = (a + 0x6d2b79f5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

export function genericBeads(name, size) {
  // A lumpy, asymmetric globule: a modest core plus 6 lobes pushed well out,
  // so the blob reads as "a protein" rather than a ball.
  const rand = mulberry32(fnv1a(name));
  const beads = [[0, 0, 0, size * 0.2]];
  for (let i = 0; i < 6; i++) {
    const v = [rand() * 2 - 1, rand() * 2 - 1, rand() * 2 - 1];
    const k = (size * (0.3 + rand() * 0.12)) / Math.max(len(v), 0.3);
    beads.push([v[0] * k, v[1] * k, v[2] * k, size * (0.12 + rand() * 0.06)]);
  }
  const radius = Math.max(...beads.map((b) => len(b.slice(0, 3)) + b[3]));
  return { beads, radius };
}

const MC_ISOLATION = 80;
const MC_SUBTRACT = 40; // higher = tighter blobs (less merging between beads)

export function makeProtein(decl, beads, style) {
  const shape = beads && beads.beads.length ? beads : genericBeads(decl.name || 'protein', decl.params.size ?? 6);
  const mat = glowMaterial(decl.params.color, style);
  const resolution = beads && beads.beads.length > 200 ? 64 : 48;
  const mc = new MarchingCubes(resolution, mat, false, false, 200000);
  mc.isolation = MC_ISOLATION;
  const maxBall = Math.max(...shape.beads.map((b) => b[3]));
  const half = (shape.radius + maxBall * 1.6) * 1.08; // world half-size of the field cube
  mc.reset();
  for (const [x, y, z, r] of shape.beads) {
    const rn = (r * 1.3) / (2 * half);            // surface radius in field units
    const strength = (MC_ISOLATION + MC_SUBTRACT) * rn * rn;
    mc.addBall(0.5 + x / (2 * half), 0.5 + y / (2 * half), 0.5 + z / (2 * half), strength, MC_SUBTRACT);
  }
  mc.update();
  mc.scale.setScalar(half);
  const group = new THREE.Group();
  group.add(mc);
  return {
    object: group,
    worldRadius: shape.radius,
    update(s) {
      applyTransform(group, s);
      setLook(mat, s.color, s.opacity, s.glow);
    },
    dispose() { mc.geometry.dispose(); mat.dispose(); },
  };
}

// ---------------------------------------------------------------------------
// Point, arrow, label

export function makePoint(decl, _beads, style) {
  const mat = glowMaterial(decl.params.color, style);
  const mesh = new THREE.Mesh(new THREE.SphereGeometry(0.8, 20, 14), mat);
  return {
    object: mesh,
    update(s) {
      applyTransform(mesh, s);
      mesh.visible = mesh.visible && !!decl.params.show;
      setLook(mat, s.color, s.opacity, s.glow);
    },
    dispose() { mesh.geometry.dispose(); mat.dispose(); },
  };
}

export function makeArrow(decl, _beads, style) {
  const mat = glowMaterial(decl.params.color, style);
  const group = new THREE.Group();
  const shaft = new THREE.Mesh(new THREE.CylinderGeometry(0.15, 0.15, 1, 10), mat);
  const head = new THREE.Mesh(new THREE.ConeGeometry(0.5, 1.2, 14), mat);
  group.add(shaft, head);
  const up = new THREE.Vector3(0, 1, 0);
  return {
    object: group,
    update(s) {
      group.visible = s.visible && s.opacity > 0.002;
      setLook(mat, s.color, s.opacity, s.glow);
      const a = new THREE.Vector3(...s.ends.from), b = new THREE.Vector3(...s.ends.to);
      const dir = b.clone().sub(a);
      const L = dir.length();
      if (L < 1e-4) { group.visible = false; return; }
      dir.normalize();
      const start = a.clone().addScaledVector(dir, L * 0.08);
      const end = b.clone().addScaledVector(dir, -L * 0.08);
      const bodyLen = Math.max(0.01, start.distanceTo(end) - 1.2);
      const q = new THREE.Quaternion().setFromUnitVectors(up, dir);
      shaft.quaternion.copy(q);
      shaft.scale.set(1, bodyLen, 1);
      shaft.position.copy(start).addScaledVector(dir, bodyLen / 2);
      head.quaternion.copy(q);
      head.position.copy(end).addScaledVector(dir, -0.6);
    },
    dispose() { shaft.geometry.dispose(); head.geometry.dispose(); mat.dispose(); },
  };
}

// Labels are drawn by the 2D overlay; they have no 3D presence.
export function makeLabel() {
  const object = new THREE.Group();
  return { object, update() {}, dispose() {} };
}

