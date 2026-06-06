"use strict";

// ---------- interaction line shades (grayscale, no hue) ----------
const COLORS = {
  hydrogen_bond: 0x222222,
  hydrophobic: 0x999999,
  salt_bridge: 0x555555,
  metal_coordination: 0x000000,
  aromatic: 0x777777,
};
const TYPE_LABEL = {
  hydrogen_bond: "H-bond",
  hydrophobic: "Hydrophobic",
  salt_bridge: "Salt bridge",
  metal_coordination: "Metal",
  aromatic: "Aromatic",
};
const WATER = ["HOH", "WAT", "DOD", "H2O", "SOL"];

// ---------- app state ----------
const state = {
  pdbId: null,
  pdbData: null,
  meta: null,
  components: [],
  selectedComp: null,
  profile: null,
  report: null,
  chains: null,
  proteinAtomCount: null,
  chemical: null,
  dockPose: null,
  dockData: null,
  screen: null,
  methods: null,
  pockets: [],
  dockSite: null,
  pocketView: null,
  evolution: null,
  colorByConservation: false,
  showCoupling: false,
  showDivergence: false,
  viewer: null,
  showSurface: false,
  showLines: true,
};

// ---------- dom helpers ----------
const $ = (sel) => document.querySelector(sel);
const el = (tag, cls, html) => {
  const e = document.createElement(tag);
  if (cls) e.className = cls;
  if (html !== undefined) e.innerHTML = html;
  return e;
};
function escapeHtml(str) {
  return String(str).replace(/[&<>"']/g, (c) => (
    { "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c]
  ));
}

function setStatus(msg, kind) {
  const bar = $("#statusBar");
  bar.className = "status-bar" + (kind ? " " + kind : "");
  bar.innerHTML = kind === "busy" ? `<span class="spinner"></span>${msg}` : msg;
}
function switchTab(name) {
  document.querySelectorAll(".tab").forEach((t) =>
    t.classList.toggle("active", t.dataset.tab === name)
  );
  document.querySelectorAll(".panel").forEach((p) =>
    p.classList.toggle("active", p.dataset.panel === name)
  );
  if (name === "viewer" && state.viewer) {
    setTimeout(() => {
      state.viewer.resize();
      state.viewer.render();
    }, 30);
  }
  if (name === "report") compileReport();
}

async function getJSON(url) {
  const resp = await fetch(url);
  const data = await resp.json();
  if (!resp.ok || data.error) throw new Error(data.error || `HTTP ${resp.status}`);
  return data;
}

// ================= structure loading =================
// One box for everything: a 4-char PDB ID loads directly; anything else
// is treated as a name and searched against RCSB.
function smartLoad(value) {
  const v = (value || "").trim();
  if (!v) {
    setStatus("Enter a PDB ID (e.g. 1HSG) or a protein name (e.g. insulin).", "error");
    return;
  }
  if (/^[0-9A-Za-z]{4}$/.test(v)) {
    loadStructure(v);
  } else {
    searchPDB(v);
  }
}

async function loadStructure(pdbId) {
  pdbId = (pdbId || "").trim().toUpperCase();
  if (!/^[0-9A-Z]{4}$/.test(pdbId)) {
    setStatus("PDB IDs are 4 characters, e.g. 1HSG.", "error");
    return;
  }
  setStatus(`Fetching ${pdbId} from RCSB and parsing atoms…`, "busy");
  $("#loadBtn").disabled = true;
  $("#searchResults").innerHTML = "";
  try {
    const data = await getJSON(`/api/analyze?pdb=${pdbId}`);
    state.pdbId = pdbId;
    state.pdbData = data.pdb_data;
    state.meta = data.metadata;
    state.components = data.components;
    state.chains = data.chains;
    state.proteinAtomCount = data.protein_atom_count;
    state.selectedComp = null;
    state.profile = null;
    state.report = null;
    state.dockData = null;
    state.screen = null;
    state.dockPose = null;
    state.methods = null;
    state.pockets = [];
    state.dockSite = null;
    state.pocketView = null;
    state.evolution = null;
    state.colorByConservation = false;
    state.showCoupling = false;
    state.showDivergence = false;
    const evoCb = $("#toggleConservation");
    if (evoCb) evoCb.checked = false;
    $("#evolutionContent").className = "empty";
    $("#evolutionContent").textContent = "No conservation analysis yet. Click “Analyze conservation”.";

    renderOverview(data);
    renderComponents(data.components);
    updateDockPocket();
    $("#pocketsContent").className = "empty";
    $("#pocketsContent").textContent = "No pockets detected yet. Click “Detect pockets”.";
    initViewer(data.pdb_data);
    switchTab("overview");

    const ligCount = data.components.filter((c) => c.kind === "ligand").length;
    setStatus(
      `Loaded ${pdbId}: ${data.protein_atom_count} protein atoms, ` +
        `${data.components.length} bound component(s) (${ligCount} ligand-like). ` +
        `Pick a molecule on the left to profile interactions.`
    );
  } catch (err) {
    setStatus(`Could not load ${pdbId}: ${err.message}`, "error");
  } finally {
    $("#loadBtn").disabled = false;
  }
}

function renderOverview(data) {
  const m = data.meta || data.metadata || {};
  const c = $("#overviewContent");
  c.className = "";
  const cell = (k, v) =>
    `<div class="meta-cell"><div class="k">${k}</div><div class="v">${
      v ?? "—"
    }</div></div>`;
  c.innerHTML = `
    <div class="title-block">
      <h3>${m.title || "Untitled structure"}</h3>
      <span class="pdbid">PDB ${m.pdb_id}</span>
    </div>
    <div class="meta-grid">
      ${cell("Method", m.experimental_method)}
      ${cell("Resolution", m.resolution_A ? m.resolution_A + " Å" : "—")}
      ${cell("Released", m.deposited ? m.deposited.slice(0, 10) : "—")}
      ${cell("Protein chains", data.chains ? data.chains.join(", ") : "—")}
      ${cell("Protein atoms", data.protein_atom_count)}
      ${cell("Bound components", data.components.length)}
      ${cell("Mol. weight", m.molecular_weight_kDa ? m.molecular_weight_kDa + " kDa" : "—")}
      ${cell("Deposited atoms", m.deposited_atom_count)}
    </div>`;
}

function renderComponents(components) {
  const card = $("#componentCard");
  const list = $("#componentList");
  list.innerHTML = "";
  if (!components.length) {
    card.hidden = false;
    list.innerHTML = `<p class="hint">No bound ligands, ions, or metals in this structure (apo). Try 1HSG or 1CA2.</p>`;
    return;
  }
  card.hidden = false;
  components.forEach((comp) => {
    const row = el("div", "comp");
    row.dataset.index = comp.index;
    row.innerHTML = `
      <span class="tag ${comp.kind}">${comp.kind}</span>
      <div>
        <div class="comp-label">${comp.res_name}</div>
        <div class="comp-sub">${comp.chain}/${comp.res_seq} · ${comp.atom_count} atoms</div>
      </div>`;
    row.addEventListener("click", () => selectComponent(comp.index));
    list.appendChild(row);
  });
}

// ================= 3D viewer =================
function initViewer(pdbData) {
  const host = $("#viewer3d");
  host.innerHTML = "";
  state.viewer = $3Dmol.createViewer(host, { backgroundColor: "#ffffff" });
  rebuildScene(true);
}

function xyz(arr) {
  return { x: arr[0], y: arr[1], z: arr[2] };
}

function drawInteractionLines(list) {
  const v = state.viewer;
  list.forEach((it) => {
    v.addCylinder({
      start: xyz(it.ligand_atom.xyz),
      end: xyz(it.protein_atom.xyz),
      radius: 0.05,
      color: COLORS[it.type] || 0x000000,
      fromCap: 1,
      toCap: 1,
    });
  });
}

function rebuildScene(resetZoom) {
  const v = state.viewer;
  if (!v || !state.pdbData) return;
  v.removeAllModels();
  v.removeAllShapes();
  v.removeAllSurfaces();
  v.removeAllLabels();

  // model 0: protein + crystallographic hetero
  v.addModel(state.pdbData, "pdb");
  v.setStyle({}, { cartoon: { color: "#b3b3b3", opacity: 0.85 } });
  if (state.colorByConservation && state.evolution && state.evolution.residues) {
    applyConservationColors(v);
  }
  v.setStyle({ hetflag: true }, { stick: { radius: 0.16, color: "#666666" } });
  v.setStyle({ resn: WATER }, {});

  const comp = state.selectedComp;
  if (comp) {
    v.setStyle(
      { chain: comp.chain, resi: comp.res_seq },
      { stick: { radius: 0.26, color: "#1a1a1a" }, sphere: { scale: 0.24, color: "#1a1a1a" } }
    );
  }
  if (state.profile && state.showLines) {
    state.profile.contact_residues.forEach((r) =>
      v.addStyle({ chain: r.chain, resi: r.res_seq }, { stick: { radius: 0.12, color: "#808080" } })
    );
    drawInteractionLines(state.profile.interactions);
  }

  // model 1: docked pose (if any)
  if (state.dockPose) {
    v.addModel(state.dockPose.pdb, "pdb");
    v.setStyle(
      { model: 1 },
      { stick: { radius: 0.3, color: "#000000" }, sphere: { scale: 0.22, color: "#000000" } }
    );
    state.dockPose.profile.contact_residues.forEach((r) =>
      v.addStyle({ model: 0, chain: r.chain, resi: r.res_seq }, { stick: { radius: 0.12, color: "#808080" } })
    );
    if (state.showLines) drawInteractionLines(state.dockPose.profile.interactions);
  }

  // detected pocket: translucent sphere at cavity center + lining residues
  if (state.pocketView) {
    const p = state.pocketView;
    const r = Math.cbrt((3 * p.volume_A3) / (4 * Math.PI));
    v.addSphere({
      center: { x: p.center[0], y: p.center[1], z: p.center[2] },
      radius: Math.max(2.0, Math.min(r, 9.0)),
      color: 0x333333,
      opacity: 0.25,
    });
    p.lining_residues.forEach((rr) =>
      v.addStyle({ chain: rr.chain, resi: rr.res_seq }, { stick: { radius: 0.12, color: "#808080" } })
    );
  }

  // coevolution network: lines between co-evolving residue pairs + hub sticks
  if (state.showCoupling && state.evolution && state.evolution.coupling_reliable) {
    state.evolution.coupling_pairs.forEach((p) => {
      if (p.xyz_i && p.xyz_j) {
        v.addCylinder({
          start: xyz(p.xyz_i),
          end: xyz(p.xyz_j),
          radius: 0.1,
          color: 0x444444,
          fromCap: 1,
          toCap: 1,
        });
      }
    });
    state.evolution.coupling_hubs.forEach((h) => {
      const resi = parseInt(String(h.res).replace(/\D/g, ""), 10);
      if (!isNaN(resi))
        v.addStyle({ chain: h.chain, resi }, { stick: { radius: 0.2, color: "#000000" } });
    });
  }

  // ancestral-divergence: highlight residues that differ from the family consensus
  if (state.showDivergence && state.evolution && state.evolution.divergent_residues) {
    state.evolution.divergent_residues.forEach((r) => {
      const resi = parseInt(String(r.res).replace(/\D/g, ""), 10);
      if (!isNaN(resi))
        v.addStyle(
          { chain: r.chain, resi },
          { stick: { radius: 0.22, color: "#1a1a1a" }, sphere: { scale: 0.22, color: "#1a1a1a" } }
        );
    });
  }

  if (state.showSurface) {
    v.addSurface($3Dmol.SurfaceType.VDW, { opacity: 0.5, color: "#d0d0d0" }, { model: 0, hetflag: false });
  }
  if (resetZoom) v.zoomTo();
  v.render();
}

function updateDockPocket() {
  const el = $("#dockPocket");
  if (state.dockSite) {
    el.textContent = state.dockSite.label;
    el.classList.add("set");
  } else {
    el.textContent = "Pick a bound molecule (step 2) or detect a pocket to set the target.";
    el.classList.remove("set");
  }
}

// ================= interaction profiling =================
async function selectComponent(index) {
  const comp = state.components.find((c) => c.index === index);
  if (!comp) return;
  document.querySelectorAll(".comp").forEach((r) =>
    r.classList.toggle("active", Number(r.dataset.index) === index)
  );
  setStatus(`Profiling atomic interactions for ${comp.label}…`, "busy");
  try {
    const data = await getJSON(
      `/api/interactions?pdb=${state.pdbId}&comp=${index}`
    );
    state.selectedComp = comp;
    state.profile = data.profile;
    state.report = data.report;
    state.dockPose = null;
    state.methods = null;
    state.pocketView = null;
    state.dockSite = { type: "comp", index: comp.index, label: comp.label };
    updateDockPocket();

    renderInteractions(data.profile);
    compileReport();
    rebuildScene();
    if (state.viewer) {
      state.viewer.zoomTo({ chain: comp.chain, resi: comp.res_seq });
      state.viewer.zoom(0.55, 800);
    }
    switchTab("interactions");
    const t = data.profile.interaction_total;
    setStatus(
      `${comp.label}: ${t} atomic interactions across ` +
        `${data.profile.contact_residue_count} residues. See the 3D viewer for contact geometry.`
    );
  } catch (err) {
    setStatus(`Interaction analysis failed: ${err.message}`, "error");
  }
}

function renderInteractions(profile) {
  const c = $("#interactionContent");
  c.className = "";
  c.innerHTML = interactionsHTML(profile);
}

function interactionsHTML(profile) {
  const counts = profile.counts;
  const chip = (type) =>
    `<div class="count-chip">
      <div class="n">${counts[type] || 0}</div>
      <div class="l">${TYPE_LABEL[type]}</div>
    </div>`;

  const rows = profile.interactions
    .map((it) => {
      const la = it.ligand_atom;
      const pa = it.protein_atom;
      return `<tr>
        <td><span class="pill ${it.type}">${TYPE_LABEL[it.type]}</span></td>
        <td>${la.name} <span style="color:var(--muted)">(${la.element})</span></td>
        <td>${pa.res_name}${pa.res_seq} · ${pa.name} <span style="color:var(--muted)">${pa.chain}</span></td>
        <td>${it.distance} Å</td>
      </tr>`;
    })
    .join("");

  const resChips = profile.contact_residues
    .map(
      (r) =>
        `<div class="res-chip"><b>${r.res_name}${r.res_seq}</b> <span class="rd">${r.chain} · ${r.total}× · ${r.min_distance}Å</span></div>`
    )
    .join("");

  return `
    <div class="counts-row">
      ${chip("hydrogen_bond")}
      ${chip("hydrophobic")}
      ${chip("salt_bridge")}
      ${chip("metal_coordination")}
      ${chip("aromatic")}
    </div>
    <div class="section-h">Atomic contacts — ${profile.component.label}</div>
    <table class="data">
      <thead><tr><th>Type</th><th>Ligand atom</th><th>Protein atom</th><th>Distance</th></tr></thead>
      <tbody>${rows || `<tr><td colspan="4">No heavy-atom contacts within cutoffs.</td></tr>`}</tbody>
    </table>
    <div class="section-h">Binding-site residues (perturbation hot spots)</div>
    <div class="res-chips">${resChips || "—"}</div>`;
}

// ================= docking =================
async function runDock() {
  const chem = $("#dockChemInput").value.trim();
  switchTab("docking");
  if (!state.pdbId) {
    setStatus("Load a structure first.", "error");
    return;
  }
  if (!chem) {
    setStatus("Enter a chemical to dock.", "error");
    return;
  }
  if (!state.dockSite) {
    setStatus("Pick a bound molecule (step 2) or a detected pocket as the target.", "error");
    return;
  }
  setStatus(
    `Docking ${chem} into ${state.dockSite.label} — Monte-Carlo search, a few seconds…`,
    "busy"
  );
  $("#dockBtn").disabled = true;
  try {
    const siteParam =
      state.dockSite.type === "comp"
        ? `comp=${state.dockSite.index}`
        : `pocket=${state.dockSite.index}`;
    const data = await getJSON(
      `/api/dock?pdb=${state.pdbId}&chem=${encodeURIComponent(chem)}&${siteParam}`
    );
    state.dockPose = { pdb: data.pose_pdb, profile: data.profile, report: data.report };
    state.dockData = data;
    state.methods = data.methods;
    renderDocking(data);
    compileReport();
    rebuildScene(false);
    if (state.viewer) {
      state.viewer.zoomTo({ model: 1 });
      state.viewer.zoom(0.5, 600);
    }
    setStatus(
      `Docked ${data.chemical.name}: score ${data.docking.score}, ` +
        `${data.profile.interaction_total} predicted interactions across ` +
        `${data.profile.contact_residue_count} residues.`
    );
  } catch (err) {
    setStatus(`Docking failed: ${err.message}`, "error");
    $("#dockingContent").className = "empty";
    $("#dockingContent").textContent = err.message;
  } finally {
    $("#dockBtn").disabled = false;
  }
}

// Reproducibility / methods block (docking-literature reporting standard).
function methodsHTML(m) {
  if (!m) return "";
  const r = m.receptor || {};
  const b = m.box || {};
  const s = m.search || {};
  const cut = m.interaction_cutoffs_A || {};
  const row = (k, v) =>
    v === undefined || v === null || v === ""
      ? ""
      : `<tr><td class="mk">${k}</td><td>${escapeHtml(String(v))}</td></tr>`;
  const lig = m.ligand
    ? row("Ligand", `PubChem CID ${m.ligand.cid} · ${m.ligand.conformer} conformer · ${m.ligand.n_heavy_atoms} heavy atoms · ${m.ligand.flexibility}`)
    : "";
  return `
    <details class="methods">
      <summary>Methods &amp; reproducibility</summary>
      <table class="methods-table">
        ${row("Tool", m.tool)}
        ${row("Run", m.run_utc)}
        ${row("Receptor", [r.pdb_id, r.title].filter(Boolean).join(" — "))}
        ${row("Experiment", [r.method, r.resolution_A ? r.resolution_A + " Å" : null].filter(Boolean).join(", "))}
        ${row("Receptor prep", m.receptor_prep)}
        ${row("Site", m.site)}
        ${row("Box", `center [${(b.center||[]).join(", ")}] · ${b.edge_A} Å edge · ${b.grid_spacing_A} Å grid · ±${b.translation_search_A} Å search`)}
        ${row("Scoring", m.scoring)}
        ${row("Search", `${s.algorithm} · ${s.seeds} seeds × ${s.mc_steps} steps · random seed ${s.random_seed}`)}
        ${lig}
        ${row("Interaction cutoffs", `H-bond ≤${cut.hydrogen_bond} · salt ≤${cut.salt_bridge} · hydrophobic ≤${cut.hydrophobic} · metal ≤${cut.metal_coordination} · aromatic ≤${cut.aromatic_centroid} Å`)}
      </table>
      <div class="disclaimer" style="margin-top:10px">${m.disclaimer || ""}</div>
    </details>`;
}

function renderDocking(data) {
  const c = $("#dockingContent");
  c.className = "";
  const d = data.docking;
  const srcNote =
    data.chemical.coord_source === "2d"
      ? " ⚠ only a 2D conformer was available — pose is approximate"
      : "";
  const box = (n, l) => `<div class="score-box"><div class="n">${n}</div><div class="l">${l}</div></div>`;
  const rmsd = d.redock_rmsd;
  const rmsdBox = rmsd != null ? box(rmsd + " Å", "Redock RMSD vs crystal") : "";
  const rmsdNote =
    rmsd != null
      ? ` Redock RMSD to the crystallographic ligand is <b>${rmsd} Å</b> (under ~2 Å = pose reproduced).`
      : "";
  c.innerHTML = `
    <div class="score-banner">
      ${box(d.score, "Docking score (lower = better)")}
      ${box(d.ligand_efficiency, "Per-atom score")}
      ${box(data.chemical.n_heavy_atoms, "Heavy atoms")}
      ${box(data.profile.contact_residue_count, "Contact residues")}
      ${rmsdBox}
    </div>
    <div class="dock-note">
      Docked <b>${data.chemical.name}</b> (CID ${data.chemical.cid}, ${data.chemical.formula || ""})
      into <b>${data.pocket.label}</b>${srcNote}. Black sticks in the 3D viewer show the predicted pose.${rmsdNote}
    </div>
    ${interactionsHTML(data.profile)}
    ${pharmacologyHTML(data.pharmacology)}
    ${methodsHTML(data.methods)}`;
}

// ================= comprehensive session report =================
// Build the report once as structured sections; render both HTML and text from
// the same data so the on-screen report and the .txt export always match.
function buildReportSections() {
  const s = [];
  const m = state.meta || {};
  const pct = (x) => (x == null ? "—" : Math.round(x * 100) + "%");

  // --- Overview ---
  if (state.pdbId) {
    s.push({
      title: "Structure overview",
      rows: [
        ["PDB", `${m.pdb_id || state.pdbId}${m.title ? " — " + m.title : ""}`],
        ["Method", [m.experimental_method, m.resolution_A ? m.resolution_A + " Å" : null].filter(Boolean).join(", ") || "—"],
        ["Released", m.deposited ? m.deposited.slice(0, 10) : "—"],
        ["Chains", state.chains ? state.chains.join(", ") : "—"],
        ["Protein atoms", state.proteinAtomCount ?? "—"],
        ["Bound components", (state.components || []).length],
      ],
      list: (state.components || []).map((c) => `${c.res_name} (${c.kind}, ${c.chain}/${c.res_seq})`),
    });
  }

  // --- Bound-ligand interaction analysis ---
  if (state.selectedComp && state.profile) {
    const p = state.profile;
    s.push({
      title: `Interaction analysis — ${p.component.label}`,
      rows: Object.entries(p.counts).map(([k, v]) => [TYPE_LABEL[k] || k, v]),
      lines: [
        `Contacts: ${p.interaction_total} across ${p.contact_residue_count} residues`,
        "Binding-site residues: " + (p.contact_residues || []).slice(0, 10).map((r) => r.res_name + r.res_seq).join(", "),
      ],
      hypotheses: state.report ? state.report.hypotheses : [],
    });
  }

  // --- Pockets ---
  if ((state.pockets || []).length) {
    const evoMap = {};
    if (state.evolution && state.evolution.pocket_conservation)
      state.evolution.pocket_conservation.forEach((pc) => (evoMap[pc.index] = pc));
    s.push({
      title: `Detected pockets (${state.pockets.length})`,
      table: {
        head: ["#", "Tier", "Volume Å³", "Druggability", "Enclosure", "Conservation", "Assessment"],
        rows: state.pockets.map((p) => {
          const e = evoMap[p.index] || {};
          return [
            "#" + (p.index + 1), p.tier, p.volume_A3, p.score, p.enclosure + "/7",
            e.mean_conservation == null ? "—" : e.mean_conservation,
            e.label || "—",
          ];
        }),
      },
    });
  }

  // --- Evolution ---
  if (state.evolution && state.evolution.available !== false) {
    const e = state.evolution;
    const rows = [
      ["Pfam family", `${e.pfam} — ${e.family_name || ""}`],
      ["Homologs", e.n_sequences],
      ["Coverage", `${e.mapped_residues}/${e.target_length} (${pct(e.coverage)})`],
      ["Consensus identity", pct(e.consensus_identity)],
      ["Coevolution confidence", `${pct(e.coupling_confidence)} (${e.coupling_reliable ? "reliable" : "suppressed — too weak"})`],
    ];
    const lines = [
      "Most conserved: " + (e.top_conserved || []).slice(0, 8).map((r) => `${r.res}(${r.conservation})`).join(", "),
    ];
    if (e.coupling_reliable && e.coupling_pairs.length)
      lines.push("Top co-evolving pairs: " + e.coupling_pairs.slice(0, 6).map((p) => `${p.res_i}–${p.res_j}(${p.distance_A}Å)`).join(", "));
    if ((e.divergent_residues || []).length)
      lines.push("Derived (vs consensus): " + e.divergent_residues.slice(0, 8).map((r) => `${r.res} ${r.from}→${r.to}`).join(", "));
    s.push({ title: "Evolutionary analysis", rows, lines });
  }

  // --- Docking ---
  if (state.dockData) {
    const d = state.dockData;
    const rows = [
      ["Chemical", `${d.chemical.name} (CID ${d.chemical.cid}, ${d.chemical.formula || ""})`],
      ["Site", d.pocket.label],
      ["Docking score", `${d.docking.score} (per-atom ${d.docking.ligand_efficiency})`],
      ["Predicted interactions", `${d.profile.interaction_total} across ${d.profile.contact_residue_count} residues`],
    ];
    if (d.docking.redock_rmsd != null) rows.push(["Redock RMSD", `${d.docking.redock_rmsd} Å vs crystal ligand`]);
    s.push({
      title: "Docking",
      rows,
      hypotheses: d.report ? d.report.hypotheses : [],
      pharmacology: d.pharmacology,
    });
  }

  // --- Chemical lookup ---
  if (state.chemical) {
    const d = state.chemical;
    const dl = d.druglikeness || {};
    const r = d.rules || {};
    const ab = r.absorption || {};
    const passmark = (x) => (!x || x.pass == null ? "—" : x.pass ? "pass" : "fail");
    s.push({
      title: `Chemical — ${d.iupac_name || d.query}`,
      rows: [
        ["CID / formula", `${d.cid} · ${d.molecular_formula || ""}`],
        ["MW / XLogP / TPSA", `${d.molecular_weight || "—"} · ${d.xlogp ?? "—"} · ${d.tpsa ?? "—"}`],
        ["Lipinski", dl.drug_like ? "drug-like" : `${dl.violation_count} violation(s)`],
        ["Veber / Egan / Lead", `${passmark(r.veber)} / ${passmark(r.egan)} / ${passmark(r.lead_like)}`],
        ["GI absorption / BBB", `${ab.gi_absorption || "—"} / ${ab.bbb_permeant || "—"}`],
      ],
      pharmacology: d.pharmacology,
    });
  }

  // --- Screening ---
  if (state.screen && (state.screen.results || []).length) {
    const ok = state.screen.results.filter((r) => !r.error);
    s.push({
      title: `Virtual screen — ${state.screen.site}`,
      table: {
        head: ["Rank", "Chemical", "Score", "Per-atom", "H-bonds", "Salt", "Contacts"],
        rows: ok.map((r) => ["#" + r.rank, r.query, r.score, r.ligand_efficiency, r.counts.hydrogen_bond, r.counts.salt_bridge, r.contact_residue_count]),
      },
    });
  }

  return s;
}

// Cross-module synthesis: observations that combine findings.
function synthesizeFindings() {
  const out = [];
  const e = state.evolution;
  const pk = state.pockets || [];

  if (pk.length && e && e.pocket_conservation) {
    const cons = e.pocket_conservation.filter((p) => /conserved \(likely functional\)/.test(p.label));
    const allo = e.pocket_conservation.filter((p) => /allosteric/.test(p.label));
    const spec = e.pocket_conservation.filter((p) => p.specificity_candidate);
    if (cons.length) out.push(`Pocket(s) ${cons.map((p) => "#" + (p.index + 1)).join(", ")} are evolutionarily conserved → likely functional/orthosteric site(s).`);
    if (allo.length) out.push(`Pocket(s) ${allo.map((p) => "#" + (p.index + 1)).join(", ")} are coupling-enriched but not conserved → candidate allosteric/cryptic control site(s).`);
    if (spec.length) out.push(`Pocket(s) ${spec.map((p) => "#" + (p.index + 1)).join(", ")} are lined by residues this protein diverged from the family consensus → possible lineage-specific binding specialization.`);
  }
  if (state.dockData) {
    const d = state.dockData;
    const ph = d.pharmacology;
    if (ph && ph.match && (ph.match.level === "uniprot" || ph.match.level === "name")) {
      const act = ph.match.best_activity ? ` (measured ${ph.match.best_activity.type} ${ph.match.best_activity.relation} ${ph.match.best_activity.value_nM} nM)` : "";
      out.push(`Docked ${d.chemical.name} is a KNOWN modulator of this target${act} — the predicted pose is consistent with experimental pharmacology.`);
    } else if (ph) {
      out.push(`Docked ${d.chemical.name} has ChEMBL pharmacology but no measured activity against this exact target — treat the predicted pose as a novel hypothesis.`);
    }
    if (d.docking.redock_rmsd != null && d.docking.redock_rmsd <= 2.5)
      out.push(`Redock RMSD ${d.docking.redock_rmsd} Å indicates the docking protocol reproduces the crystallographic pose for this system.`);
  }
  if (!out.length) out.push("Run more modules (pockets + evolution + docking) to generate cross-analysis synthesis.");
  return out;
}

function compileReport() {
  const c = $("#reportContent");
  if (!state.pdbId) {
    c.className = "empty";
    c.textContent = "Load a structure and run some analyses, then compile a report.";
    return;
  }
  const sections = buildReportSections();
  const synth = synthesizeFindings();

  const rowsHTML = (rows) =>
    `<table class="methods-table">${rows.map(([k, v]) => `<tr><td class="mk">${k}</td><td>${v}</td></tr>`).join("")}</table>`;
  const tableHTML = (t) =>
    `<table class="data"><thead><tr>${t.head.map((h) => `<th>${h}</th>`).join("")}</tr></thead><tbody>${t.rows.map((r) => `<tr>${r.map((x) => `<td>${x}</td>`).join("")}</tr>`).join("")}</tbody></table>`;

  let html = `<div class="report-summary">Comprehensive analysis report for <b>${state.meta.pdb_id || state.pdbId}</b>${state.meta.title ? " — " + state.meta.title : ""}. Generated by AtomScope from the analyses run this session.</div>`;
  html += `<div class="section-h">Integrated synthesis</div><ul class="hyp-list">${synth.map((x) => `<li>${x}</li>`).join("")}</ul>`;

  sections.forEach((sec) => {
    html += `<div class="section-h">${sec.title}</div>`;
    if (sec.rows) html += rowsHTML(sec.rows);
    if (sec.list && sec.list.length) html += `<div class="hint">${sec.list.join(" · ")}</div>`;
    if (sec.lines) html += sec.lines.map((l) => `<div class="hint" style="margin-top:4px">${l}</div>`).join("");
    if (sec.table) html += tableHTML(sec.table);
    if (sec.pharmacology) html += pharmacologyHTML(sec.pharmacology);
    if (sec.hypotheses && sec.hypotheses.length)
      html += `<ul class="hyp-list" style="margin-top:8px">${sec.hypotheses.map((h) => `<li>${h}</li>`).join("")}</ul>`;
  });

  if (state.methods) html += methodsHTML(state.methods);
  html += `<div class="disclaimer">Research-only. AtomScope interactions, docking, pockets, and coevolution are geometric/empirical/statistical heuristics from a single static structure and family alignment — not affinities, structures, or clinical guidance. Validate with orthogonal evidence.</div>`;

  c.className = "";
  c.innerHTML = html;
}

// ================= batch virtual screen =================
async function runScreen() {
  const raw = $("#screenInput").value.trim();
  if (!state.pdbId) {
    setStatus("Load a structure first.", "error");
    return;
  }
  if (!raw) {
    setStatus("Enter one or more chemicals to screen.", "error");
    return;
  }
  if (!state.dockSite) {
    setStatus("Set a docking target first (pick a molecule or a pocket).", "error");
    return;
  }
  setStatus(
    `Screening into ${state.dockSite.label} — docking each chemical on a shared grid, ~10–20 s…`,
    "busy"
  );
  $("#screenBtn").disabled = true;
  try {
    const siteParam =
      state.dockSite.type === "comp"
        ? `comp=${state.dockSite.index}`
        : `pocket=${state.dockSite.index}`;
    const data = await getJSON(
      `/api/screen?pdb=${state.pdbId}&chems=${encodeURIComponent(raw)}&${siteParam}`
    );
    state.screen = data;
    renderScreen(data);
    const ok = data.results.filter((r) => !r.error).length;
    setStatus(`Screened ${data.results.length} chemical(s) into ${data.site}; ${ok} docked and ranked by fit.`);
  } catch (err) {
    setStatus(`Screen failed: ${err.message}`, "error");
  } finally {
    $("#screenBtn").disabled = false;
  }
}

function renderScreen(data) {
  const c = $("#screenContent");
  const rows = data.results
    .map((r) => {
      if (r.error) {
        return `<tr class="bad-row"><td>—</td><td>${r.query}</td><td colspan="6">${r.error}</td><td></td></tr>`;
      }
      return `<tr>
        <td class="rank-cell">#${r.rank}</td>
        <td>${r.query}</td>
        <td>${r.formula || ""}</td>
        <td>${r.score}</td>
        <td>${r.ligand_efficiency}</td>
        <td>${r.counts.hydrogen_bond}</td>
        <td>${r.counts.salt_bridge}</td>
        <td>${r.contact_residue_count}</td>
        <td><button class="mini view-btn" data-chem="${r.query}">View pose</button></td>
      </tr>`;
    })
    .join("");
  c.innerHTML = `
    <div class="section-h" style="margin-top:18px">Ranking — ${data.site}</div>
    <table class="data">
      <thead><tr><th>Rank</th><th>Chemical</th><th>Formula</th><th>Score</th><th>Per-atom</th><th>H-bonds</th><th>Salt</th><th>Contacts</th><th></th></tr></thead>
      <tbody>${rows}</tbody>
    </table>
    <div class="hint">Lower score = better predicted fit. “View pose” docks that chemical singly and shows the pose in the 3D viewer.</div>
    ${methodsHTML(data.methods)}`;
  c.querySelectorAll(".view-btn").forEach((b) =>
    b.addEventListener("click", () => {
      $("#dockChemInput").value = b.dataset.chem;
      runDock();
    })
  );
}

// ================= pocket detection =================
async function detectPockets() {
  if (!state.pdbId) {
    setStatus("Load a structure first.", "error");
    return;
  }
  switchTab("pockets");
  setStatus("Detecting pockets (LIGSITE geometric cavity scan)… a few seconds.", "busy");
  $("#detectBtn").disabled = true;
  try {
    const data = await getJSON(`/api/pockets?pdb=${state.pdbId}`);
    state.pockets = data.pockets;
    renderPockets(data.pockets);
    setStatus(
      `Found ${data.count} candidate pocket(s). Click “Show in 3D” to view, or “Dock here” to dock a chemical into a cavity.`
    );
  } catch (err) {
    setStatus(`Pocket detection failed: ${err.message}`, "error");
  } finally {
    $("#detectBtn").disabled = false;
  }
}

function pocketCard(p) {
  const res = p.lining_residues
    .slice(0, 12)
    .map((r) => r.res_name + r.res_seq)
    .join(", ");
  const more = p.lining_residue_count > 12 ? ", …" : "";
  const ss = p.subscores || {};
  const pct = (x) => Math.round((x || 0) * 100) + "%";
  const tier = p.tier || "pocket";
  return `<div class="pocket-card" id="pk-card-${p.index}">
    <div class="pocket-top">
      <span class="pocket-rank">#${p.index + 1}</span>
      <span class="tier-badge ${tier}" title="pocket = geometric cavity · ligandable = big &amp; enclosed enough · druggable = also chemically favourable">${tier}</span>
      <div class="scorebar"><i style="width:${p.score}%"></i></div>
      <div class="pocket-metrics">
        <span>druggability <b>${p.score}</b></span>
        <span>volume <b>${p.volume_A3} Å³</b></span>
        <span>enclosure <b>${p.enclosure}/7</b></span>
        <span>lining <b>${p.lining_residue_count}</b></span>
      </div>
      <div class="pocket-actions">
        <button class="mini" id="pk-show-${p.index}">Show in 3D</button>
        <button class="mini" id="pk-dock-${p.index}">Dock here</button>
      </div>
    </div>
    <div class="pocket-sub">hydrophobic <b>${pct(ss.hydrophobicity)}</b> · polar/charged <b>${pct(ss.polarity)}</b> · aromatic <b>${pct(ss.aromaticity)}</b></div>
    <div class="pocket-res">Lining residues: ${res}${more}</div>
  </div>`;
}

function renderPockets(pk) {
  const c = $("#pocketsContent");
  if (!pk.length) {
    c.className = "empty";
    c.textContent = "No enclosed cavities found (the surface may be open/flat).";
    return;
  }
  c.className = "";
  const tiers = { druggable: 0, ligandable: 0, pocket: 0 };
  pk.forEach((p) => (tiers[p.tier || "pocket"] += 1));
  const legend = `<div class="tier-legend">
    <span class="tier-badge druggable">druggable</span> ${tiers.druggable}
    &nbsp;·&nbsp; <span class="tier-badge ligandable">ligandable</span> ${tiers.ligandable}
    &nbsp;·&nbsp; <span class="tier-badge pocket">pocket</span> ${tiers.pocket}
    <span class="hint" style="display:block;margin-top:6px">Tiers (heuristic): <b>pocket</b> = geometric cavity · <b>ligandable</b> = large &amp; enclosed enough for a small molecule · <b>druggable</b> = also chemically favourable (hydrophobic, enclosed, not too polar).</span>
  </div>`;
  c.innerHTML = legend + `<div class="pocket-list">${pk.map(pocketCard).join("")}</div>`;
  pk.forEach((p) => {
    c.querySelector(`#pk-show-${p.index}`).addEventListener("click", () => showPocket(p.index));
    c.querySelector(`#pk-dock-${p.index}`).addEventListener("click", () => dockIntoPocket(p.index));
  });
}

function showPocket(index) {
  const p = state.pockets.find((x) => x.index === index);
  if (!p) return;
  state.pocketView = p;
  document.querySelectorAll(".pocket-card").forEach((el) => el.classList.remove("active"));
  const card = $(`#pk-card-${index}`);
  if (card) card.classList.add("active");
  rebuildScene(false);
  if (state.viewer) {
    const resi = p.lining_residues.map((r) => r.res_seq);
    state.viewer.zoomTo({ resi });
    state.viewer.zoom(0.85, 500);
  }
  switchTab("viewer");
}

function dockIntoPocket(index) {
  const p = state.pockets.find((x) => x.index === index);
  if (!p) return;
  state.pocketView = p;
  state.dockSite = {
    type: "pocket",
    index,
    label: `detected pocket #${index + 1} (${p.volume_A3} Å³)`,
  };
  updateDockPocket();
  switchTab("docking");
  setStatus(`Target set to pocket #${index + 1}. Enter a chemical and click Dock.`);
  $("#dockChemInput").focus();
}

// ================= evolution / conservation =================
async function runEvolution() {
  if (!state.pdbId) {
    setStatus("Load a structure first.", "error");
    return;
  }
  switchTab("evolution");
  setStatus("Fetching the Pfam family alignment, scoring conservation and computing coevolution — up to ~30s…", "busy");
  $("#evoBtn").disabled = true;
  try {
    const data = await getJSON(`/api/evolution?pdb=${state.pdbId}`);
    if (!data.available) {
      state.evolution = null;
      $("#evolutionContent").className = "empty";
      $("#evolutionContent").textContent = data.reason || "No conservation data available.";
      setStatus("No Pfam family / alignment found for this structure.", "error");
      return;
    }
    state.evolution = data;
    renderEvolution(data);
    setStatus(
      `Conservation from ${data.n_sequences} Pfam homologs (${data.pfam}); ` +
        `${Math.round(data.coverage * 100)}% of residues mapped.`
    );
  } catch (err) {
    setStatus(`Conservation analysis failed: ${err.message}`, "error");
  } finally {
    $("#evoBtn").disabled = false;
  }
}

function cellEvo(k, v) {
  return `<div class="meta-cell"><div class="k">${k}</div><div class="v">${v}</div></div>`;
}

function renderEvolution(d) {
  const c = $("#evolutionContent");
  c.className = "";
  const top = d.top_conserved
    .map((r) => `<span class="res-chip"><b>${r.res}</b> <span class="rd">${r.conservation}</span></span>`)
    .join("");
  const pk = (d.pocket_conservation || [])
    .map(
      (p) => `<tr>
        <td>#${p.index + 1}</td>
        <td>${p.tier || "—"}</td>
        <td>${p.volume_A3} Å³</td>
        <td>${p.mean_conservation == null ? "—" : p.mean_conservation}</td>
        <td>${p.divergent_lining || 0}${p.specificity_candidate ? " ★" : ""}</td>
        <td>${p.label}</td>
      </tr>`
    )
    .join("");
  c.innerHTML = `
    <div class="meta-grid">
      ${cellEvo("Pfam family", d.pfam + " — " + (d.family_name || ""))}
      ${cellEvo("Homologs", d.n_sequences)}
      ${cellEvo("Residues mapped", `${d.mapped_residues} / ${d.target_length} (${Math.round(d.coverage * 100)}%)`)}
      ${cellEvo("UniProt", d.uniprot || "—")}
    </div>
    <button id="evoColorBtn" class="primary" style="margin:16px 0">Color structure by conservation →</button>
    <div class="section-h">Most conserved residues</div>
    <div class="res-chips">${top}</div>
    ${coevolutionHTML(d)}
    ${divergenceHTML(d)}
    <div class="section-h">Pocket conservation &amp; divergence</div>
    <table class="data">
      <thead><tr><th>Pocket</th><th>Tier</th><th>Volume</th><th>Mean conservation</th><th>Divergent lining</th><th>Assessment</th></tr></thead>
      <tbody>${pk || `<tr><td colspan="6">No pockets.</td></tr>`}</tbody>
    </table>
    <div class="hint">★ = specificity candidate: pocket lined by residues this protein has diverged from the conserved family consensus (possible lineage-specific binding specialization).</div>
    <div class="disclaimer">Conservation is a family-MSA signal (Shannon entropy across ${d.n_sequences} Pfam homologs, 0 = variable, 1 = invariant) — not phylogenetic ancestral reconstruction. Conserved pocket-lining residues suggest functional importance.</div>`;
  $("#evoColorBtn").addEventListener("click", () => {
    state.colorByConservation = true;
    const cb = $("#toggleConservation");
    if (cb) cb.checked = true;
    rebuildScene(false);
    switchTab("viewer");
  });
  const netBtn = $("#netBtn");
  if (netBtn)
    netBtn.addEventListener("click", () => {
      state.showCoupling = true;
      rebuildScene(false);
      switchTab("viewer");
    });
  const divBtn = $("#divBtn");
  if (divBtn)
    divBtn.addEventListener("click", () => {
      state.showDivergence = true;
      rebuildScene(false);
      switchTab("viewer");
    });
}

function divergenceHTML(d) {
  const idPct = Math.round((d.consensus_identity || 0) * 100);
  const rows = (d.divergent_residues || [])
    .map(
      (r) => `<tr>
        <td>${r.res} <span style="color:var(--muted)">${r.chain}</span></td>
        <td>${r.from} → ${r.to}</td>
        <td>${r.conservation}</td>
      </tr>`
    )
    .join("");
  return `
    <div class="section-h">Ancestral divergence (family consensus)</div>
    <div class="hint" style="margin-top:0">This protein is <b>${idPct}%</b> identical to its family consensus (an ancestral-like "average" sequence). The positions below are where it has <b>diverged at otherwise-conserved columns</b> — lineage-specific "derived" substitutions that often determine specialized function or specificity. This is an ASR-style proxy, not a phylogenetic ancestral reconstruction.</div>
    ${rows ? `<button id="divBtn" class="primary" style="margin:12px 0">Highlight divergent residues in 3D →</button>` : ""}
    <table class="data" style="margin-top:8px">
      <thead><tr><th>Residue</th><th>Consensus → this protein</th><th>Conservation</th></tr></thead>
      <tbody>${rows || `<tr><td colspan="3">No notable divergence at conserved positions.</td></tr>`}</tbody>
    </table>`;
}

function coevolutionHTML(d) {
  const conf = Math.round((d.coupling_confidence || 0) * 100);
  if (!d.coupling_reliable) {
    return `
      <div class="section-h">Evolutionary coupling (MIp)</div>
      <div class="pharm-match miss">Coevolution signal is too weak to trust for this family — only ${conf}% of the top co-evolving pairs are spatial contacts in this structure (a deeper / more divergent alignment is needed). Network and allosteric flags are suppressed rather than shown as noise.</div>`;
  }
  const hubs = d.coupling_hubs
    .map((h) => `<span class="res-chip"><b>${h.res}</b> <span class="rd">deg ${h.degree}</span></span>`)
    .join("");
  const pairs = d.coupling_pairs
    .slice(0, 15)
    .map(
      (p) => `<tr>
        <td>${p.res_i} – ${p.res_j}</td>
        <td>${p.mip}</td>
        <td>${p.distance_A == null ? "—" : p.distance_A + " Å"}</td>
      </tr>`
    )
    .join("");
  return `
    <div class="section-h">Evolutionary coupling (MIp) — confidence ${conf}%</div>
    <div class="hint" style="margin-top:0">${conf}% of top co-evolving pairs are spatial contacts here, so the signal is trusted. Co-evolving residue networks frequently mark functional or allosteric couplings (the EVcouplings/Gremlin idea, approximated by APC-corrected mutual information).</div>
    <button id="netBtn" class="primary" style="margin:12px 0">Show coevolution network in 3D →</button>
    <div class="section-h">Coupling hubs</div>
    <div class="res-chips">${hubs || "—"}</div>
    <table class="data" style="margin-top:10px">
      <thead><tr><th>Co-evolving pair</th><th>MIp</th><th>CA–CA distance</th></tr></thead>
      <tbody>${pairs}</tbody>
    </table>`;
}

function applyConservationColors(v) {
  const groups = {};
  state.evolution.residues.forEach((r) => {
    const c = r.conservation;
    let shade;
    if (c == null) shade = "#e2e2e2";
    else if (c >= 0.7) shade = "#1a1a1a";
    else if (c >= 0.5) shade = "#555555";
    else if (c >= 0.3) shade = "#999999";
    else shade = "#cccccc";
    const key = r.chain + "|" + shade;
    (groups[key] = groups[key] || { chain: r.chain, shade, resi: [] }).resi.push(r.res_seq);
  });
  Object.values(groups).forEach((g) =>
    v.setStyle({ chain: g.chain, resi: g.resi }, { cartoon: { color: g.shade } })
  );
}

// ================= chemical lookup =================
async function lookupChemical(q) {
  q = (q || "").trim();
  if (!q) return;
  const ctx = state.pdbId
    ? ` and cross-referencing ChEMBL vs ${state.pdbId}…`
    : "…";
  setStatus(`Looking up "${q}" in PubChem${ctx}`, "busy");
  switchTab("chemical");
  try {
    const pdbParam = state.pdbId ? `&pdb=${state.pdbId}` : "";
    const data = await getJSON(`/api/chemical?q=${encodeURIComponent(q)}${pdbParam}`);
    state.chemical = data;
    if (!$("#dockChemInput").value.trim()) $("#dockChemInput").value = q;
    renderChemical(data);
    setStatus(`Loaded chemical: ${data.molecular_formula || q} (CID ${data.cid}).`);
  } catch (err) {
    setStatus(`Chemical lookup failed: ${err.message}`, "error");
    $("#chemicalContent").className = "empty";
    $("#chemicalContent").textContent = err.message;
  }
}

function renderChemical(d) {
  const c = $("#chemicalContent");
  c.className = "";
  const dl = d.druglikeness || {};
  const dlClass = dl.drug_like ? "ok" : "warn";
  const dlBadge = dl.drug_like
    ? `<span class="badge-ok">drug-like ✓</span>`
    : `<span class="badge-warn">${dl.violation_count} rule violation(s)</span>`;
  const violations = (dl.violations || []).length
    ? `<div class="hint">${dl.violations.join(" · ")}</div>`
    : "";

  let chembl = "";
  if (d.chembl) {
    chembl = `<div class="druglike ${d.chembl.max_phase === 4 ? "ok" : ""}">
      <div><b>ChEMBL:</b> ${d.chembl.pref_name || d.chembl.chembl_id} —
      <span class="${d.chembl.max_phase === 4 ? "badge-ok" : ""}">${d.chembl.development_status}</span></div>
      ${d.chembl.url ? `<div class="hint"><a class="ext" href="${d.chembl.url}" target="_blank">View in ChEMBL ↗</a></div>` : ""}
    </div>`;
  }

  const prop = (k, v) =>
    `<div class="meta-cell"><div class="k">${k}</div><div class="v">${v ?? "—"}</div></div>`;

  c.innerHTML = `
    <div class="chem-head">
      ${d.image_url ? `<img class="chem-img" src="${d.image_url}" alt="2D structure" />` : ""}
      <div class="chem-info">
        <h3>${d.iupac_name || d.query}</h3>
        <div class="pdbid" style="font-family:monospace">CID ${d.cid} · ${d.molecular_formula || ""}</div>
        ${d.pubchem_url ? `<div class="hint"><a class="ext" href="${d.pubchem_url}" target="_blank">View in PubChem ↗</a></div>` : ""}
        <div class="druglike ${dlClass}">
          <div>${dl.rule || "Druglikeness"}: ${dlBadge}</div>
          ${violations}
        </div>
        ${chembl}
      </div>
    </div>
    ${pharmacologyHTML(d.pharmacology)}
    <div class="section-h">Druglikeness rules</div>
    ${rulesHTML(d.rules)}
    <div class="section-h">Physicochemical properties</div>
    <div class="meta-grid">
      ${propCell("Molecular weight", d.molecular_weight ? d.molecular_weight + " g/mol" : "—")}
      ${propCell("XLogP", d.xlogp)}
      ${propCell("TPSA", d.tpsa ? d.tpsa + " Å²" : "—")}
      ${propCell("H-bond donors", d.h_bond_donors)}
      ${propCell("H-bond acceptors", d.h_bond_acceptors)}
      ${propCell("Rotatable bonds", d.rotatable_bonds)}
      ${propCell("Formal charge", d.formal_charge)}
    </div>
    <div class="section-h">SMILES</div>
    <div class="smiles">${d.smiles || "—"}</div>`;
}

function propCell(k, v) {
  return `<div class="meta-cell"><div class="k">${k}</div><div class="v">${v ?? "—"}</div></div>`;
}

// Druglikeness rule chips (Veber / Egan / lead-like + BOILED-Egg absorption).
function rulesHTML(rules) {
  if (!rules) return `<div class="hint">—</div>`;
  const mark = (r) =>
    !r || r.pass === null ? "—" : r.pass ? "✓" : "✗";
  const cls = (r) => (!r || r.pass === null ? "" : r.pass ? "pass" : "fail");
  const ab = rules.absorption || {};
  return `<div class="rule-row">
    <span class="rule-chip ${cls(rules.veber)}" title="${rules.veber ? rules.veber.criteria : ""}">Veber ${mark(rules.veber)}</span>
    <span class="rule-chip ${cls(rules.egan)}" title="${rules.egan ? rules.egan.criteria : ""}">Egan ${mark(rules.egan)}</span>
    <span class="rule-chip ${cls(rules.lead_like)}" title="${rules.lead_like ? rules.lead_like.criteria : ""}">Lead-like ${mark(rules.lead_like)}</span>
    <span class="rule-chip" title="${ab.model || ""}">GI absorption: <b>${ab.gi_absorption || "—"}</b></span>
    <span class="rule-chip" title="${ab.model || ""}">BBB: <b>${ab.bbb_permeant || "—"}</b></span>
  </div>`;
}

// Curated ChEMBL pharmacology + measured activity vs the loaded protein.
function pharmacologyHTML(ph) {
  if (!ph) return "";
  const m = ph.match || {};
  const fmtAct = (a) =>
    a ? `${a.type} ${a.relation} ${a.value_nM} nM` : null;
  let banner = "";
  if (m.level === "uniprot" || m.level === "name") {
    const conf = m.level === "uniprot" ? "confirmed (UniProt)" : "likely (name match)";
    const act = fmtAct(m.best_activity);
    banner = `<div class="pharm-match hit">
      ✓ Known modulator of the loaded target — <b>${m.target_name}</b> <span class="conf">[${conf}]</span>
      ${act ? `· best measured <b>${act}</b>` : ""}
    </div>`;
  } else {
    banner = `<div class="pharm-match miss">No measured ChEMBL activity links this chemical to the loaded protein (it may still bind — or simply be unstudied here).</div>`;
  }
  const mechs = (ph.mechanisms || [])
    .map((mm) => {
      const act = fmtAct(mm.best_activity);
      return `<tr>
        <td>${mm.moa || "—"}</td>
        <td>${mm.target_name || "—"}</td>
        <td>${act || "—"}</td>
      </tr>`;
    })
    .join("");
  const mechTable = mechs
    ? `<table class="data" style="margin-top:8px">
        <thead><tr><th>Mechanism of action</th><th>Target</th><th>Best measured potency</th></tr></thead>
        <tbody>${mechs}</tbody>
       </table>`
    : `<div class="hint">No curated mechanism-of-action entries in ChEMBL.</div>`;
  return `
    <div class="section-h">Known pharmacology (ChEMBL) — ${ph.pref_name || ph.chembl_id}, ${ph.development_status}</div>
    ${banner}
    ${mechTable}
    ${ph.url ? `<div class="hint"><a class="ext" href="${ph.url}" target="_blank">View ${ph.chembl_id} in ChEMBL ↗</a></div>` : ""}`;
}

// ================= search =================
async function searchPDB(q) {
  q = (q || "").trim();
  if (!q) return;
  const box = $("#searchResults");
  box.innerHTML = `<div class="hint"><span class="spinner"></span>Searching…</div>`;
  try {
    const data = await getJSON(`/api/search?q=${encodeURIComponent(q)}`);
    if (!data.results.length) {
      box.innerHTML = `<div class="hint">No matches.</div>`;
      return;
    }
    box.innerHTML = `<div class="hint" style="margin:0 0 6px">${data.results.length} match(es) — click one to load:</div>`;
    data.results.forEach((r) => {
      const title = r.title ? escapeHtml(r.title) : "(no title)";
      const org = r.organism ? escapeHtml(r.organism) : "";
      const row = el(
        "div",
        "sr",
        `<div class="sr-top"><span class="sr-id">${r.pdb_id}</span><span class="sr-org">${org}</span></div>
         <div class="sr-title">${title}</div>`
      );
      row.addEventListener("click", () => {
        $("#pdbInput").value = r.pdb_id;
        loadStructure(r.pdb_id);
      });
      box.appendChild(row);
    });
  } catch (err) {
    box.innerHTML = `<div class="hint">Search failed: ${err.message}</div>`;
  }
}

// ================= export =================
function exportReport() {
  if (!state.pdbId) return;
  const m = state.meta || {};
  const L = [];
  const rule = (ch) => ch.repeat(60);
  L.push("AtomScope — Comprehensive Analysis Report (research-only)");
  L.push(rule("="));
  L.push(`PDB: ${m.pdb_id || state.pdbId}  ${m.title || ""}`);
  L.push(`Generated: ${new Date().toISOString().slice(0, 19).replace("T", " ")} (local)`);
  L.push("");

  L.push("INTEGRATED SYNTHESIS");
  synthesizeFindings().forEach((x, i) => L.push(`  ${i + 1}. ${x}`));
  L.push("");

  buildReportSections().forEach((sec) => {
    L.push(sec.title.toUpperCase());
    L.push(rule("-"));
    (sec.rows || []).forEach(([k, v]) => L.push(`  ${k}: ${v}`));
    (sec.list || []).forEach((x) => L.push(`  - ${x}`));
    (sec.lines || []).forEach((x) => L.push(`  ${x}`));
    if (sec.table) {
      L.push("  " + sec.table.head.join(" | "));
      sec.table.rows.forEach((r) => L.push("  " + r.join(" | ")));
    }
    if (sec.pharmacology) {
      const ph = sec.pharmacology;
      L.push(`  ChEMBL: ${ph.pref_name || ph.chembl_id} — ${ph.development_status}`);
      if (ph.match && ph.match.level !== "none")
        L.push(`  Target match (${ph.match.level}): ${ph.match.target_name}` + (ph.match.best_activity ? ` — best ${ph.match.best_activity.type} ${ph.match.best_activity.relation} ${ph.match.best_activity.value_nM} nM` : ""));
      (ph.mechanisms || []).forEach((mm) =>
        L.push(`  MoA: ${mm.moa} | ${mm.target_name}` + (mm.best_activity ? ` (${mm.best_activity.type} ${mm.best_activity.value_nM} nM)` : ""))
      );
    }
    (sec.hypotheses || []).forEach((h, i) => L.push(`  Hypothesis ${i + 1}: ${h}`));
    L.push("");
  });

  const mm = state.methods;
  if (mm) {
    L.push("METHODS & REPRODUCIBILITY");
    L.push(rule("-"));
    L.push(`  Tool: ${mm.tool}   Run: ${mm.run_utc}`);
    if (mm.receptor) L.push(`  Receptor: ${[mm.receptor.pdb_id, mm.receptor.title].filter(Boolean).join(" — ")}`);
    L.push(`  Receptor prep: ${mm.receptor_prep}`);
    L.push(`  Site: ${mm.site}`);
    if (mm.box) L.push(`  Box: center [${(mm.box.center || []).join(", ")}], ${mm.box.edge_A} A edge, ${mm.box.grid_spacing_A} A grid, +/-${mm.box.translation_search_A} A search`);
    L.push(`  Scoring: ${mm.scoring}`);
    if (mm.search) L.push(`  Search: ${mm.search.algorithm}, ${mm.search.seeds} seeds x ${mm.search.mc_steps} steps, random seed ${mm.search.random_seed}`);
    if (mm.ligand) L.push(`  Ligand: PubChem CID ${mm.ligand.cid}, ${mm.ligand.conformer} conformer, ${mm.ligand.n_heavy_atoms} heavy atoms, ${mm.ligand.flexibility}`);
    L.push("");
  }

  L.push(rule("="));
  L.push("Research-only. Heuristic predictions from a single static structure + family alignment. Not affinities or clinical guidance; validate with orthogonal evidence.");

  const blob = new Blob([L.join("\n")], { type: "text/plain" });
  const a = document.createElement("a");
  a.href = URL.createObjectURL(blob);
  a.download = `atomscope_${m.pdb_id || state.pdbId}_report.txt`;
  a.click();
  URL.revokeObjectURL(a.href);
}

// ================= wiring =================
function init() {
  $("#loadBtn").addEventListener("click", () => smartLoad($("#pdbInput").value));
  $("#pdbInput").addEventListener("keydown", (e) => {
    if (e.key === "Enter") smartLoad($("#pdbInput").value);
  });
  document.querySelectorAll("a.ex").forEach((a) =>
    a.addEventListener("click", (e) => {
      e.preventDefault();
      $("#pdbInput").value = a.dataset.pdb;
      loadStructure(a.dataset.pdb);
    })
  );
  $("#chemBtn").addEventListener("click", () => lookupChemical($("#chemInput").value));
  $("#chemInput").addEventListener("keydown", (e) => {
    if (e.key === "Enter") lookupChemical($("#chemInput").value);
  });
  $("#detectBtn").addEventListener("click", detectPockets);
  $("#evoBtn").addEventListener("click", runEvolution);
  $("#toggleConservation").addEventListener("change", (e) => {
    if (e.target.checked && !state.evolution) {
      e.target.checked = false;
      runEvolution();
      return;
    }
    state.colorByConservation = e.target.checked;
    rebuildScene(false);
  });
  $("#dockBtn").addEventListener("click", runDock);
  $("#dockChemInput").addEventListener("keydown", (e) => {
    if (e.key === "Enter") runDock();
  });
  $("#screenBtn").addEventListener("click", runScreen);
  document.querySelectorAll(".tab").forEach((t) =>
    t.addEventListener("click", () => switchTab(t.dataset.tab))
  );
  $("#toggleSurface").addEventListener("change", (e) => {
    state.showSurface = e.target.checked;
    rebuildScene();
  });
  $("#toggleLines").addEventListener("change", (e) => {
    state.showLines = e.target.checked;
    rebuildScene();
  });
  $("#resetView").addEventListener("click", () => {
    if (state.viewer) {
      state.viewer.zoomTo();
      state.viewer.render();
    }
  });
  $("#exportBtn").addEventListener("click", exportReport);
  $("#compileBtn").addEventListener("click", compileReport);
}

document.addEventListener("DOMContentLoaded", init);
