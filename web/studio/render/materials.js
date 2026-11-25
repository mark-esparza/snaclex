// SnaCleX Studio — the soft "glowing molecule" material.
//
// A standard lit material plus a Fresnel rim term: surfaces facing away from
// the viewer pick up an emissive halo in the object's own colour. On the dark
// style this gives the classic look of molecules floating in darkness; on the
// light style the rim is toned down so shapes stay crisp on a pale page.

import * as THREE from '../../vendor/three.module.min.js';

const RIM = { dark: 1.6, light: 0.45 };

export function glowMaterial(color, style = 'dark') {
  const m = new THREE.MeshStandardMaterial({
    color, roughness: 0.5, metalness: 0.05,
    emissive: color, emissiveIntensity: style === 'dark' ? 0.35 : 0.06,
    transparent: true, opacity: 1,
  });
  const uniforms = {
    rimStrength: { value: 0.4 * RIM[style] },
    rimColor: { value: new THREE.Color(color) },
  };
  m.userData.uniforms = uniforms;
  m.userData.style = style;
  m.onBeforeCompile = (shader) => {
    Object.assign(shader.uniforms, uniforms);
    shader.fragmentShader = shader.fragmentShader
      .replace('#include <common>', '#include <common>\nuniform float rimStrength;\nuniform vec3 rimColor;')
      .replace(
        '#include <emissivemap_fragment>',
        `#include <emissivemap_fragment>
        float rimF = pow(1.0 - abs(dot(normalize(vViewPosition), normal)), 2.2);
        totalEmissiveRadiance += rimColor * rimF * rimStrength;`,
      );
  };
  return m;
}

// Update colour, opacity and glow (0..1) on a material from glowMaterial().
export function setLook(m, color, opacity, glow) {
  const u = m.userData.uniforms;
  if (m.userData.color !== color) {
    m.color.set(color);
    m.emissive.set(color);
    if (u) u.rimColor.value.set(color);
    m.userData.color = color;
  }
  m.opacity = opacity;
  m.depthWrite = opacity > 0.98;
  if (u) u.rimStrength.value = glow * RIM[m.userData.style || 'dark'];
}
