// Shared test fixtures for the Studio's Node tests.

// One fixed-column PDB ATOM/HETATM record.
export function pdbLine(rec, serial, name, res, chain, seq, x, y, z, alt = ' ') {
  const atom = name.length < 4 ? ` ${name}`.padEnd(4) : name;
  return `${rec.padEnd(6)}${String(serial).padStart(5)} ${atom}${alt}${res.padStart(3)} ${chain}${String(seq).padStart(4)}    `
    + `${x.toFixed(3).padStart(8)}${y.toFixed(3).padStart(8)}${z.toFixed(3).padStart(8)}  1.00  0.00           C`;
}

// A small chain of n CA atoms spaced 3.8 Å apart along x.
export function chainPdb(chain, n, start = 0) {
  return Array.from({ length: n }, (_, i) => pdbLine('ATOM', i + 1, 'CA', 'ALA', chain, i + 1, start + i * 3.8, 0, 0)).join('\n');
}
