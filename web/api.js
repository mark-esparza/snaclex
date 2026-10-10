// Renders GET /api/docs into the API reference page.
// External, not inline: the CSP is script-src 'self' with no
// 'unsafe-inline', so an inline <script> here is silently refused by
// the browser and the page never leaves its "Loading…" state.
function esc(s) {
  return String(s == null ? "" : s).replace(/[&<>"']/g, (c) => (
    { "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c]
  ));
}
function kv(obj) {
  if (!obj || typeof obj !== "object") return esc(obj);
  return Object.entries(obj)
    .map(([k, v]) => `<div><code>${esc(k)}</code> — ${esc(v)}</div>`)
    .join("");
}
fetch("/api/docs").then((r) => r.json()).then((d) => {
  const limits = Object.entries(d.limits || {})
    .map(([k, v]) => `<tr><td class="mk">${esc(k)}</td><td>${esc(v)}</td></tr>`)
    .join("");
  const rows = (d.endpoints || []).map((e) => `
    <div class="api-ep">
      <div class="api-ep-head"><span class="api-method">${esc(e.method)}</span>
        <code>${esc(e.path)}</code></div>
      ${e.params ? `<div class="api-sub"><b>Params</b>${kv(e.params)}</div>` : ""}
      ${e.body ? `<div class="api-sub"><b>Body</b><div>${kv(e.body)}</div></div>` : ""}
      <div class="api-sub"><b>Returns</b> ${esc(e.returns)}</div>
      ${e.notes ? `<div class="api-sub muted">${esc(e.notes)}</div>` : ""}
    </div>`).join("");
  document.getElementById("apiContent").className = "";
  document.getElementById("apiContent").innerHTML = `
    <p>${esc(d.tool)} v${esc(d.version)} · base <code>${esc(d.base_url)}</code></p>
    <h2>Limits</h2>
    <table class="methods-table">${limits}</table>
    <h2>Errors</h2>
    <p>${esc(d.errors)}</p>
    <h2>Endpoints</h2>
    ${rows}`;
}).catch((err) => {
  document.getElementById("apiContent").textContent =
    "Could not load API docs: " + err.message;
});
