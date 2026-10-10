// Renders the reference benchmark cases into the help page.
//
// External, not inline: the CSP is script-src 'self' with no 'unsafe-inline',
// so an inline <script> is silently refused by the browser.
//
// The list is fetched from /api/benchmark/cases rather than hardcoded here, so
// the cases documented on this page and the cases the project benchmarks
// against cannot drift apart. This page documents them; the workbench itself
// offers no pre-loaded structures or chemicals.

function escapeHtml(str) {
  return String(str == null ? "" : str).replace(/[&<>"']/g, (c) => (
    { "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c]
  ));
}

function renderBenchCases(cases) {
  const host = document.getElementById("benchCases");
  if (!host) return;
  if (!cases || !cases.length) {
    host.textContent = "No reference cases are configured on this deployment.";
    return;
  }
  const rows = cases
    .map(
      (c) => `<tr>
        <td><code>${escapeHtml(c.pdb)}</code></td>
        <td>${escapeHtml(c.name || c.ligand || "")}</td>
        <td><code>${escapeHtml(c.ligand || "")}</code></td>
        <td>${escapeHtml(c.site || "")}</td>
      </tr>`
    )
    .join("");
  host.className = "tablebox";
  host.innerHTML = `
    <table class="data">
      <thead><tr>
        <th scope="col">PDB ID</th><th scope="col">Ligand</th>
        <th scope="col">Residue</th><th scope="col">Site</th>
      </tr></thead>
      <tbody>${rows}</tbody>
    </table>`;
}

fetch("/api/benchmark/cases")
  .then((r) => r.json())
  .then((d) => renderBenchCases(d.cases))
  .catch((err) => {
    const host = document.getElementById("benchCases");
    if (host) {
      host.textContent = "Could not load reference cases: " + err.message;
    }
  });
