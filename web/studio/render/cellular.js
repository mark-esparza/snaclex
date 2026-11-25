// SnaCleX Studio — views for cell-scale actors: cohesin-style ring,
// chromosome (one or two sister chromatids with a kinetochore), spindle fibre,
// and the cell itself (membrane, nucleus, cleavage furrow, division).
//
// These are stylised shapes, not structures: the scene's "not to scale" note
// covers them. Geometry follows the conventions in ../anchors.js so that
// anchors like chr.kinetochore land on the drawn kinetochore.

import * as THREE from '../../vendor/three.module.min.js';
import { glowMaterial, setLook } from './materials.js';
import { applyTransform } from './molecular.js';
import { KINETOCHORE_OFFSET } from '../anchors.js';

const changed = (a, b) => a === undefined || Math.abs(a - b) > 1e-3;
const Z = new THREE.Vector3(0, 0, 1);
const Y = new THREE.Vector3(0, 1, 0);

// ---------------------------------------------------------------------------
// Ring (e.g. cohesin)

const RING_TUBE = 0.35;
const MAX_GAP = Math.PI / 2; // fully open = a 90° gap
const RING_TILT = 0.5;        // radians; a ring seen exactly edge-on reads as a bar

export function makeRing(decl, _beads, style) {
  const group = new THREE.Group();   // placed and oriented
  const inner = new THREE.Group();   // holds the arc pieces
  group.add(inner);
  const mat = glowMaterial(decl.params.color, style);
  const R = decl.params.radius ?? 3;
  let last = {};
  const pieces = [];

  function rebuild(open, cleave) {
    for (const m of pieces) { inner.remove(m); m.geometry.dispose(); }
    pieces.length = 0;
    const gap = MAX_GAP * open;
    const arc = Math.PI * 2 - gap;
    // The gap is centred on +y; the ring starts just past it.
    const start = Math.PI / 2 + gap / 2;
    if (cleave < 1e-3) {
      const m = new THREE.Mesh(new THREE.TorusGeometry(R, RING_TUBE, 12, 64, arc), mat);
      m.rotation.z = start;
      pieces.push(m);
    } else {
      // Cut opposite the gap; the two halves drift apart sideways.
      const half = arc / 2 - 0.15 * cleave;
      for (const side of [0, 1]) {
        const m = new THREE.Mesh(new THREE.TorusGeometry(R, RING_TUBE, 12, 32, half), mat);
        m.rotation.z = side === 0 ? start : start + arc / 2 + 0.15 * cleave;
        const dir = side === 0 ? -1 : 1;
        m.position.set(dir * 3 * cleave, -1.2 * cleave, 0);
        pieces.push(m);
      }
    }
    for (const m of pieces) inner.add(m);
  }

  return {
    object: group,
    update(s) {
      applyTransform(group, s);
      if (s.ends && s.ends.axis) {
        // Encircle whatever it is around: ring axis (local z) along that axis.
        group.quaternion.setFromUnitVectors(Z, new THREE.Vector3(...s.ends.axis));
        group.rotateX(RING_TILT);
      }
      const open = s.p.open ?? 0, cleave = s.p.cleave ?? 0;
      if (changed(last.open, open) || changed(last.cleave, cleave)) {
        rebuild(open, cleave);
        last = { open, cleave };
      }
      setLook(mat, s.color, s.opacity * (1 - 0.5 * cleave), s.glow);
    },
    dispose() { for (const m of pieces) m.geometry.dispose(); mat.dispose(); },
  };
}

// ---------------------------------------------------------------------------
// Chromosome

const SISTER_OFFSET = 0.7;

// Lathe profile of one chromatid: a rounded rod pinched at the centromere.
function chromatidGeometry(length, radius) {
  const pts = [];
  const n = 48;
  for (let i = 0; i <= n; i++) {
    const y = -length / 2 + (length * i) / n;
    const e = Math.abs(y) / (length / 2);
    const cap = Math.sqrt(Math.max(0, 1 - e ** 8));
    const waist = 1 - 0.45 * Math.exp(-((y / 0.9) ** 2));
    pts.push(new THREE.Vector2(Math.max(0.001, radius * cap * waist), y));
  }
  return new THREE.LatheGeometry(pts, 28);
}

export function makeChromosome(decl, _beads, style) {
  const group = new THREE.Group();
  const mat = glowMaterial(decl.params.color, style);
  const kMat = glowMaterial('#f2f2f2', style);
  const sisters = Math.round(decl.params.sisters ?? 2);
  const meshes = [];
  for (let i = 0; i < sisters; i++) {
    const m = new THREE.Mesh(undefined, mat);
    meshes.push(m);
    group.add(m);
  }
  const kinetochore = new THREE.Mesh(new THREE.CylinderGeometry(0.55, 0.55, 0.18, 24), kMat);
  kinetochore.rotation.z = Math.PI / 2; // disc faces local -x
  kinetochore.position.set(-KINETOCHORE_OFFSET, 0, 0);
  group.add(kinetochore);
  let last = {};

  return {
    object: group,
    update(s) {
      applyTransform(group, s);
      setLook(mat, s.color, s.opacity, s.glow);
      setLook(kMat, '#f2f2f2', s.opacity, Math.max(s.glow, 0.6));
      const condense = s.p.condense ?? 1;
      const sep = (s.p.separate ?? 0) * (s.p.gap ?? 6);
      if (changed(last.condense, condense)) {
        const L = (decl.params.length ?? 12) * (1.4 - 0.4 * condense);
        const r = 1.0 * (0.6 + 0.4 * condense);
        const geo = chromatidGeometry(L, r);
        for (const m of meshes) { if (m.geometry) m.geometry.dispose(); m.geometry = geo; }
        last.condense = condense;
      }
      if (sisters === 2) {
        meshes[0].position.set(-SISTER_OFFSET - sep / 2, 0, 0);
        meshes[1].position.set(SISTER_OFFSET + sep / 2, 0, 0);
      }
      kinetochore.position.set(-KINETOCHORE_OFFSET - sep / 2, 0, 0); // matches anchors.js
    },
    dispose() {
      if (meshes[0] && meshes[0].geometry) meshes[0].geometry.dispose();
      kinetochore.geometry.dispose(); mat.dispose(); kMat.dispose();
    },
  };
}

// ---------------------------------------------------------------------------
// Spindle fibre (a microtubule bundle drawn as one glowing filament)

export function makeFibre(decl, _beads, style) {
  const mat = glowMaterial(decl.params.color, style);
  const mesh = new THREE.Mesh(new THREE.CylinderGeometry(0.12, 0.12, 1, 8), mat);
  return {
    object: mesh,
    update(s) {
      setLook(mat, s.color, s.opacity, Math.max(s.glow, 0.5));
      const a = new THREE.Vector3(...s.ends.from), b = new THREE.Vector3(...s.ends.tip);
      const dir = b.clone().sub(a);
      const L = dir.length();
      mesh.visible = s.visible && s.opacity > 0.002 && L > 1e-3;
      if (!mesh.visible) return;
      mesh.quaternion.setFromUnitVectors(Y, dir.normalize());
      mesh.scale.set(s.scale, L, s.scale);
      mesh.position.copy(a).add(b).multiplyScalar(0.5);
    },
    dispose() { mesh.geometry.dispose(); mat.dispose(); },
  };
}

// ---------------------------------------------------------------------------
// Cell

const MEMBRANE_OPACITY = 0.08;
const NUCLEUS_OPACITY = 0.3;

export function makeCell(decl, _beads, style) {
  const R = decl.params.radius ?? 40;
  const group = new THREE.Group();
  const mat = glowMaterial(decl.params.color, style);
  mat.side = THREE.DoubleSide;
  const nucMat = glowMaterial('#7f8cff', style);
  const base = new THREE.SphereGeometry(R, 64, 40);
  const basePos = base.attributes.position.array.slice();
  const membrane = new THREE.Mesh(base, mat);
  const nucleus = new THREE.Mesh(new THREE.SphereGeometry(R * 0.45, 40, 28), nucMat);
  // Daughter cells after division (each with a little less than the volume).
  const daughterGeo = new THREE.SphereGeometry(R * 0.8, 48, 32);
  const daughters = [new THREE.Mesh(daughterGeo, mat), new THREE.Mesh(daughterGeo, mat)];
  group.add(membrane, nucleus, ...daughters);
  const hasNucleus = decl.params.nucleus !== false;
  let lastPinch;

  function pinchTo(pinch) {
    // Cleavage furrow in the plane x = 0: squeeze y and z near the equator.
    const pos = base.attributes.position.array;
    for (let i = 0; i < pos.length; i += 3) {
      const x = basePos[i];
      const k = 1 - 0.8 * pinch * Math.exp(-((x / (0.28 * R)) ** 2));
      pos[i] = x;
      pos[i + 1] = basePos[i + 1] * k;
      pos[i + 2] = basePos[i + 2] * k;
    }
    base.attributes.position.needsUpdate = true;
    base.computeVertexNormals();
  }

  return {
    object: group,
    update(s) {
      applyTransform(group, s);
      const pinch = s.p.pinch ?? 0, split = s.p.split ?? 0;
      if (changed(lastPinch, pinch)) { pinchTo(pinch); lastPinch = pinch; }
      setLook(mat, s.color, s.opacity * MEMBRANE_OPACITY, Math.max(s.glow, 0.6));
      setLook(nucMat, '#7f8cff', s.opacity * NUCLEUS_OPACITY * (s.p.envelope ?? 1), s.glow);
      mat.depthWrite = false;
      nucMat.depthWrite = false;
      nucleus.visible = hasNucleus && (s.p.envelope ?? 1) > 0.01 && split < 1e-3;
      const divided = split > 1e-3;
      membrane.visible = !divided;
      const off = R * 0.5 + ((s.p.gap ?? 50) * split) / 2;
      daughters[0].visible = daughters[1].visible = divided;
      daughters[0].position.set(-off, 0, 0);
      daughters[1].position.set(off, 0, 0);
    },
    dispose() { base.dispose(); daughterGeo.dispose(); nucleus.geometry.dispose(); mat.dispose(); nucMat.dispose(); },
  };
}
