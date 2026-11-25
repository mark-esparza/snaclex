// SnaCleX Studio — the Export control. Filled in by the video-export task;
// until then the control explains that export is not available yet.
export function initExport() {
  const sel = document.getElementById('export');
  sel.replaceChildren(new Option('Not available yet', ''));
  sel.disabled = true;
}
