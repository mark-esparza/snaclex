// SnaCleX Studio — the 3D stage: one view per actor, a camera, lights, and
// the dark or light scene style. Stage.apply(frameState) poses everything for
// one moment; Stage.render() draws it. The stage never decides *what* happens
// — that is the timeline's job — it only shows a FrameState.

import * as THREE from '../../vendor/three.module.min.js';
import { makeHelix, makeStrand, makeProtein, makePoint, makeArrow, makeLabel } from './molecular.js';

export const VIEW_BUILDERS = {
  helix: makeHelix,
  strand: makeStrand,
  protein: makeProtein,
  point: makePoint,
  arrow: makeArrow,
  label: makeLabel,
};

export function registerViews(map) { Object.assign(VIEW_BUILDERS, map); }

const STYLES = {
  // fog: exponential density × (1 / camera distance), so depth cueing looks
  // the same whether the camera is framing a protein or a whole cell.
  dark: { background: '#05060a', fog: 0.42, hemi: ['#cdd8ff', '#2a1830', 1.6], key: 2.2 },
  light: { background: '#f4f1ea', fog: 0, hemi: ['#ffffff', '#b8b2a5', 1.4], key: 1.6 },
};

export class Stage {
  constructor(canvas, { style = 'dark' } = {}) {
    this.canvas = canvas;
    this.renderer = new THREE.WebGLRenderer({ canvas, antialias: true, preserveDrawingBuffer: true });
    this.renderer.outputColorSpace = THREE.SRGBColorSpace;
    this.scene = new THREE.Scene();
    this.camera = new THREE.PerspectiveCamera(35, 16 / 9, 0.1, 8000);
    this.hemi = new THREE.HemisphereLight('#ffffff', '#000000', 1);
    this.key = new THREE.DirectionalLight('#ffffff', 1.5);
    this.key.position.set(0.6, 1, 0.8);
    this.scene.add(this.hemi, this.key, this.key.target);
    this.views = new Map();
    this.width = canvas.width || 1280;
    this.height = canvas.height || 720;
    this.setStyle(style);
  }

  setStyle(style) {
    this.style = STYLES[style] ? style : 'dark';
    const s = STYLES[this.style];
    this.scene.background = new THREE.Color(s.background);
    this.scene.fog = s.fog ? new THREE.FogExp2(s.background, s.fog) : null;
    this.hemi.color.set(s.hemi[0]);
    this.hemi.groundColor.set(s.hemi[1]);
    this.hemi.intensity = s.hemi[2];
    this.key.intensity = s.key;
  }

  build(tl, beads = new Map()) {
    this.clear();
    this.tl = tl;
    this.setStyle(tl.style);
    for (const name of tl.order) {
      const decl = { name, ...tl.actors[name] };
      const make = VIEW_BUILDERS[decl.type];
      if (!make) continue;
      const view = make(decl, beads.get(name) || null, this.style);
      this.views.set(name, view);
      this.scene.add(view.object);
    }
  }

  clear() {
    for (const v of this.views.values()) {
      this.scene.remove(v.object);
      v.dispose();
    }
    this.views.clear();
  }

  apply(f) {
    for (const [name, view] of this.views) {
      const s = f.actors[name];
      if (s) view.update(s, f, this.tl);
    }
    const c = f.camera;
    this.camera.fov = c.fov;
    this.camera.position.set(c.pos[0], c.pos[1], c.pos[2]);
    this.camera.up.set(0, 1, 0);
    this.camera.lookAt(c.target[0], c.target[1], c.target[2]);
    this.camera.updateProjectionMatrix();
    this.camera.updateMatrixWorld();
    if (this.scene.fog) {
      const d = Math.hypot(c.pos[0] - c.target[0], c.pos[1] - c.target[1], c.pos[2] - c.target[2]);
      this.scene.fog.density = STYLES[this.style].fog / Math.max(d, 1);
    }
    // Key light rides with the camera, from above-left, aimed at the target,
    // so whatever the camera looks at stays lit.
    const dist = Math.max(1, this.camera.position.distanceTo(this.key.target.position.set(c.target[0], c.target[1], c.target[2])));
    this.key.position.copy(this.camera.position).add(new THREE.Vector3(-0.5, 1, 0.2).multiplyScalar(dist * 0.6));
  }

  render() { this.renderer.render(this.scene, this.camera); }

  setSize(width, height, pixelRatio = 1) {
    this.width = width;
    this.height = height;
    this.renderer.setPixelRatio(pixelRatio);
    this.renderer.setSize(width, height, false);
    this.camera.aspect = width / height;
    this.camera.updateProjectionMatrix();
  }

  // World point -> canvas pixels (CSS size), or null when behind the camera.
  project(v) {
    const p = new THREE.Vector3(v[0], v[1], v[2]).project(this.camera);
    if (p.z > 1 || p.z < -1) return null;
    return [((p.x + 1) / 2) * this.width, ((1 - p.y) / 2) * this.height];
  }

  // Radius used for label offsets: real-structure proteins know their size.
  viewRadius(name) {
    const v = this.views.get(name);
    return v && v.worldRadius ? v.worldRadius : null;
  }

  dispose() {
    this.clear();
    this.renderer.dispose();
  }
}
