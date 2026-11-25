// SnaCleX Studio — the Export control: lists only the formats this browser
// can encode, runs the export with a progress bar and a Cancel button, and
// downloads the finished file.

import { supportedFormats, exportVideo, formatInfo, EXPORT_SIZES } from './record.js';

const LABEL = { webm: 'WebM', mp4: 'MP4' };
const NO_WEBCODECS = 'Video export needs a browser with WebCodecs (current Chrome, Edge, Firefox or Safari).';

export async function initExport({ stage, state, slug, download, setBusy }) {
  const sel = document.getElementById('export');
  const status = document.getElementById('exportStatus');
  const progress = document.getElementById('exportProgress');
  const cancel = document.getElementById('exportCancel');

  const formats = await supportedFormats();
  sel.replaceChildren(new Option(formats.length ? 'Choose…' : 'Unavailable', ''));
  for (const format of formats) {
    for (const size of Object.keys(EXPORT_SIZES)) {
      sel.add(new Option(`${size} ${LABEL[format]}`, `${format}:${size}`));
    }
  }
  sel.disabled = !formats.length;
  if (!formats.length) status.textContent = NO_WEBCODECS;

  let controller = null;
  cancel.addEventListener('click', () => { if (controller) controller.abort(); });

  sel.addEventListener('change', async () => {
    const value = sel.value;
    sel.value = '';
    if (!value) return;
    if (!state.tl) {
      status.textContent = 'Nothing to export yet: write or load a scene first.';
      return;
    }
    const [format, size] = value.split(':');
    const [width, height] = EXPORT_SIZES[size];
    const tl = state.tl;
    controller = new AbortController();
    setBusy(true);
    sel.disabled = true;
    cancel.hidden = false;
    progress.hidden = false;
    progress.value = 0;
    status.textContent = 'Exporting… 0%';
    try {
      const blob = await exportVideo({
        timeline: tl, stage, width, height, format,
        signal: controller.signal,
        onProgress: (p) => { progress.value = p; status.textContent = `Exporting… ${Math.round(p * 100)}%`; },
      });
      const name = `${slug(tl.title)}-${size}.${formatInfo(format).ext}`;
      download(blob, name);
      status.textContent = `Saved ${name} (${(blob.size / 1e6).toFixed(1)} MB).`;
    } catch (e) {
      status.textContent = e && e.name === 'AbortError' ? 'Export cancelled.' : `Export failed: ${e && e.message ? e.message : e}`;
    } finally {
      controller = null;
      cancel.hidden = true;
      progress.hidden = true;
      sel.disabled = false;
      setBusy(false);
    }
  });
}
