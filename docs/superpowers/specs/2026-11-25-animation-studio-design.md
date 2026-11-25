# SnaCleX Animation Studio: design

Status: draft for approval
Date: 2026-11-25

## Intent and the why

The user wants a **reusable tool** for making short educational animations of
molecular and cellular processes. The look they have in mind is the classic
"molecules on a dark background" style of biology-education animation: a
glowing DNA helix, soft blob-like proteins drifting in and acting on it, and a
camera that moves through the scene.

The first animation to make with it is **meiotic nondisjunction leading to
trisomy 21**. It shows cohesin holding chromatids together, separase cutting
cohesin, spindle fibres attaching to kinetochores, a Meiosis I error, and the
n+1 gamete that becomes a 2n+1 zygote. The tool must not be built *for* that
one story, though. It must be equally usable for transcription, replication,
translation, mitosis and similar topics.

SnaCleX already loads real structures from the PDB. That connection is what
sets this tool apart from a generic animation app. Protein shapes can come from
real structures, and every video can credit the structures it used.

## Chosen approach: scene scripts

Three approaches were considered:

1. **Scene script (chosen).** The user writes, or asks Claude to draft, a short
   plain-text script. The script declares what is on screen (actors) and what
   happens when (actions on a timeline). The studio shows a live preview next
   to the script and records the result to a video file.
   - This is the fastest route to real videos and the most flexible.
   - It fits a "describe it in English, Claude drafts the script" workflow.
   - Cost: the user edits text rather than dragging things around.
2. **Visual drag-and-keyframe editor.** This is the most intuitive option, but
   it is roughly 3–4 times the work, and most of that goes into the editing UI
   rather than the animation. Deferred. It can be layered on top of the script
   format later, because the editor would just write scripts.
3. **Fixed templates with knobs** (pick "transcription" and adjust speed and
   colour). Quick to use, but it only covers topics someone has already built,
   which contradicts "a tool I can reuse". Rejected.

The user picked approach 1.

## Goals

- A new page in SnaCleX, `/studio.html`, with a script editor on the left, a 3D
  preview on the right, and a timeline scrubber with play, pause and seek
  underneath.
- A small, readable **scene-script language** that a non-programmer can edit.
  Errors are reported by line number with a plain-English message.
- A library of reusable **actors** and **actions**, listed below, that covers
  molecular-scale and cell-scale storytelling.
- The default visual style is dark background, soft glowing materials and
  blob-like proteins. It is original, not a copy of any specific existing
  video.
- **Proteins from real structures.** `protein from=7NJ1 chain=A` fetches the
  chain through SnaCleX's existing `/api/analyze` endpoint and turns it into a
  smooth coarse blob. Proteins can also be generic shapes when no structure is
  given.
- **Deterministic playback.** Scene state is a pure function of time `t`, so
  scrubbing, replaying and recording always produce the same frames.
- **Video export** to a file the user can upload to YouTube or a slide deck.
  WebM is the baseline; MP4 is offered where the browser can encode it. Export
  sizes are 720p and 1080p.
- **Automatic credits.** An end card or caption can be generated listing every
  PDB entry used, as in "Structures: 7NJ1 (separase–securin), 6WG3 (cohesin)…".
  This keeps SnaCleX's provenance ethos.
- **Save and load scripts** as `.txt` files on the user's computer. Ship
  bundled examples, starting with:
  - `nondisjunction-meiosis-1.txt`, the trisomy 21 story
  - `transcription-basic.txt`, a second topic that proves the tool is general

## Non-goals (for this version)

- A drag-and-drop or keyframe visual editor.
- Physics or molecular-dynamics simulation. Motion is scripted, not simulated.
- Narration or audio. Users add voice-over in a normal video editor.
- Photoreal or ray-traced rendering, and server-side rendering.
- User accounts or storing scripts on the server.

## Anti-goals (a failure even if it works)

- **Recreating a specific existing animation**, such as reproducing DNALC's
  videos shot for shot or copying their look exactly. The style genre is fine;
  copies are not.
- **Breaking SnaCleX's zero-dependency, no-build-step deployment**, or loosening
  the Content-Security-Policy to load scripts from new CDNs.
- **Implying a stylised animation is a measurement.** The tool must make it
  easy to label scenes "not to scale". Real structures are credited, and
  generic shapes are never presented as real structures.
- **A script language that only a programmer can write.**

## Constraints found in the codebase

- SnaCleX is Python stdlib only (`requirements.txt` is empty by design). The
  server serves `web/` as static files (`server.py`, `_serve_static`).
- The CSP (`server.py`, `_CSP`) allows scripts only from `'self'`,
  `https://3Dmol.org` and `blob:`. It also has `connect-src 'self'`.
  - The 3D engine must therefore be **vendored** into `web/vendor/` (a pinned
    `three.module.min.js` with its MIT licence added to `NOTICE`), not loaded
    from a CDN.
  - Video recording must not need any external host.
- 3Dmol.js is good for viewing structures, but it is not suited to custom
  animated geometry. Three.js (vendored) is used for the studio. 3Dmol stays
  as it is for the existing workbench.
- `/api/analyze?pdb=ID&chain=X` returns `pdb_data`. Large cryo-EM entries (for
  example 6WG3) return `too_large` unless a chain is given. Studio protein
  actors therefore always name a chain.
- CI runs `node --check web/*.js` and `python -m unittest`. E2E uses
  Playwright and Chromium against the real server with a pre-seeded offline
  cache (`e2e/test_smoke.py`). New studio JS must pass `node --check`. The
  script parser must be a separate module with no DOM or WebGL dependency, so
  it can be unit-tested in Node.

## The scene-script language (shape, not final grammar)

The script is line-based. Blank lines and lines starting with `#` are ignored.
Times are in seconds. A sketch:

```
title "Nondisjunction in Meiosis I"
size 1920x1080
duration 45s
style dark
note "Stylised; not to scale"

# --- actors ---
actor dna      = helix turns=40 color=pink
actor sep      = protein from=7NJ1 chain=A color=gold
actor sec      = protein from=7NJ1 chain=B color=orange
actor ring     = ring color=teal around=dna at=0.5
actor chr21a   = chromosome color=violet label="Chr 21 (maternal)"
actor fibre1   = fibre from=pole-left to=chr21a.kinetochore

# --- timeline ---
0s      camera wide
0-4s    sep enter from=left to=dna@0.45
4-6s    sec leave to=right
6s      caption "Securin released: separase is now active" for=3s
6-8s    ring open
8-12s   camera zoom to=ring
14-20s  fibre1 pull chr21a to=pole-left
```

Rules the builder must keep:

- Each actor line has the form `actor <name> = <type> key=value…`.
- Each timeline line has the form `<time or range> <target> <verb> key=value…`.
- Unknown words produce an error naming the line and suggesting the closest
  valid word.
- Positions can be absolute, or relative to another actor (`dna@0.45` means
  45% along the helix; `chr21a.kinetochore` names a part).
- Overlapping actions on the same actor are allowed when they touch different
  properties. Otherwise the later action wins, with a warning.
- Every action has an easing; the default is smooth (ease-in-out).

## Building blocks

Actors:

| Type | What it is | Key abilities |
|---|---|---|
| `helix` | DNA double helix along a straight or curved path | unwind or open a bubble over a range, rewind, cut, colour segments |
| `strand` | Single strand (RNA or DNA) | grow from a point, follow a path, peel away |
| `protein` | Blob from a real PDB chain, or a generic lumpy shape | enter, leave, move, ride along a helix, bind or unbind, change colour or glow |
| `ring` | Ring complex such as cohesin, optionally around a helix or chromatid pair | open, close, cleave (break apart) |
| `chromosome` | Condensed chromosome as one chromatid or sister pair (X shape), with centromere and kinetochore parts | condense, pair with homolog, separate sisters |
| `fibre` | Spindle microtubule from a pole to a target | grow, attach, pull, shrink |
| `cell` | Cell membrane, optionally with a nucleus | pinch or cleave into two, fade the nuclear envelope |
| `label`, `arrow`, `caption` | On-screen text and pointers | appear or disappear, follow an actor |

Global verbs are `enter`, `leave`, `move`, `follow`, `rotate`, `scale`, `fade`,
`color`, `glow`, `attach`, `detach` and `split`. The camera verbs are `wide`,
`zoom to`, `orbit`, `follow` and `cut to`.

## Screen layout

```
+--------------------------------------------------------------+
| SnaCleX Studio   [Examples v] [Open] [Save] [Export video v] |
+---------------------------+----------------------------------+
| 1 title "Nondisjunction…" |                                  |
| 2 duration 45s            |        3D preview (dark)         |
| 3 actor dna = helix …     |                                  |
| …                         |                                  |
| ! line 14: unknown verb   |                                  |
|   "pul" — did you mean    |                                  |
|   "pull"?                 |                                  |
+---------------------------+----------------------------------+
| > ||  0:12 ----------o----------------------------- 0:45     |
+--------------------------------------------------------------+
```

The page follows the existing SnaCleX visual language (`web/style.css`) for
its chrome. The preview canvas uses the studio's own dark style.

## Decided

- Script-based authoring (user's choice).
- A reusable tool, not a one-off video (user's choice).

## Builder's calls (proposed here, user may correct)

- It lives inside SnaCleX as `/studio.html`, so it can pull real structures
  and reuse the server. It is linked from the main navigation.
- Three.js is vendored, pinned, and loaded as an ES module from `'self'`.
- The script format is the custom line-based language above rather than JSON
  or YAML, because it is easier for a non-programmer to read and edit.
- Real-structure proteins are rendered as smooth coarse blobs. The builder
  picks the technique, such as residue-level spheres merged into a surface or
  metaballs, provided a 1,000-residue chain stays interactive in the preview.
- Recording renders frames on a fixed time step, so the video is smooth even
  on a slow machine. The builder picks the mechanism: `MediaRecorder` with
  manual frame requests, or WebCodecs where available.
- The first release ships the two example scripts named above. The
  nondisjunction script shows a **Meiosis I** error (homologs failing to
  separate) side by side with normal segregation. A Meiosis II example can
  follow. The script is checked against standard textbook descriptions.
- Parser and timeline unit tests run with Node's built-in test runner
  (`node --test`), added as one step in `ci.yml`. This adds no dependency.
- A video can be up to 5 minutes long. That is ample for short teaching clips
  and keeps memory bounded during recording.

## Left to the builder

- The exact grammar, error wording, and file layout under `web/studio/`.
- Shader and material details within the dark, soft-glow style.
- How helix unwinding and ring cleavage are animated geometrically.
- Test split: parser and timeline unit tests in Node or by an equivalent route
  that CI can run, plus a Playwright smoke test that loads an example, plays
  it, and asserts no JS errors.
