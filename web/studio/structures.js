// SnaCleX Studio — fetch real protein chains for protein actors.
//
// A protein actor with from=<PDB ID> chain=<C>[,<C>…] is fetched through
// SnaCleX's own /api/analyze endpoint (same origin; the server handles RCSB,
// caching and size limits), one request per chain. The Cα traces of all the
// chains are combined into one set of beads. If anything goes wrong the actor
// keeps a generic shape and the problem is reported as a warning on its line,
// so playback and export still work offline.

import { caCoords, coarseBeads } from './pdbcoarse.js';

const cache = new Map(); // "ID:chains" -> Promise<{beads|null, reason|null}>

export function clearStructureCache() { cache.clear(); }

async function fetchChain(id, chain, fetchFn) {
  const res = await fetchFn(`/api/analyze?pdb=${encodeURIComponent(id)}&chain=${encodeURIComponent(chain)}`);
  let body = null;
  try { body = await res.json(); } catch { body = null; }
  if (!res.ok) throw new Error((body && body.error) || `HTTP ${res.status}`);
  if (!body || body.too_large) throw new Error('structure too large');
  if (body.error) throw new Error(body.error);
  if (typeof body.pdb_data !== 'string') throw new Error('no structure data');
  const coords = caCoords(body.pdb_data, chain);
  if (!coords.length) throw new Error(`no CA atoms in chain ${chain}`);
  return coords;
}

function load(id, chains, fetchFn) {
  const key = `${id}:${chains.join(',')}`;
  if (!cache.has(key)) {
    cache.set(key, (async () => {
      try {
        const parts = [];
        for (const c of chains) parts.push(...(await fetchChain(id, c, fetchFn)));
        return { beads: coarseBeads(parts), reason: null };
      } catch (e) {
        cache.delete(key); // let a later edit retry (e.g. after going back online)
        return { beads: null, reason: e && e.message ? e.message : String(e) };
      }
    })());
  }
  return cache.get(key);
}

export async function loadStructures(scene, fetchFn = (u) => fetch(u)) {
  const beads = new Map();
  const warnings = [];
  const jobs = scene.actors
    .filter((a) => a.type === 'protein' && a.params.from)
    .map(async (a) => {
      const chains = String(a.params.chain || 'A').split(',').map((c) => c.trim()).filter(Boolean);
      const r = await load(a.params.from, chains, fetchFn);
      beads.set(a.name, r.beads);
      if (!r.beads) {
        warnings.push({
          line: a.line,
          message: `structure ${a.params.from} chain ${chains.join(',')} could not be loaded (${r.reason}); showing a generic shape`,
        });
      }
    });
  await Promise.all(jobs);
  warnings.sort((x, y) => x.line - y.line);
  return { beads, warnings };
}
