# SnaCleX Animation Studio Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `/studio.html` to SnaCleX. On that page a user writes a plain-text
scene script, watches a live dark-style 3D preview, scrubs a timeline, and
exports a WebM or MP4 video. Protein shapes can come from real PDB chains.

**Architecture:** The pure modules are `parse.js`, `schema.js`, `timeline.js`,
`anchors.js` and `pdbcoarse.js`. They turn script text into a timeline and
compute an immutable `FrameState` for any time `t`. They have no DOM and no
WebGL, and Node unit-tests them.

A thin rendering layer (`render/*.js`, Three.js) builds one view per actor and
applies a `FrameState` to it. A 2D overlay draws captions, labels and credits.
`record.js` renders frames at a fixed step and encodes them with WebCodecs plus
a vendored muxer. Everything is static files under `web/`. The server only
gains a `.txt` content type.

**Tech Stack:**
- Vanilla ES modules, with no build step.
- Three.js r170 (`three.module.min.js`) and its `MarchingCubes` add-on, vendored.
- `webm-muxer` 5.1.4 and `mp4-muxer` 5.2.2, vendored.
- WebCodecs `VideoEncoder`.
- Node 22 built-in test runner.
- Python stdlib `unittest`.
- Playwright for E2E.

**Spec:** `docs/superpowers/specs/2026-11-25-animation-studio-design.md`

## Global Constraints

- No new runtime or build dependency. `requirements.txt` stays as it is, and
  there is no npm install at run time. Vendored files are committed under
  `web/vendor/`.
- The CSP in `server.py` (`_CSP`) must not be loosened. All scripts load from
  `'self'`. No inline `<script>` and no import maps (the CSP blocks them), so
  vendored files that `import … from 'three'` are patched to
  `from './three.module.min.js'`.
- Vendored versions are pinned exactly: `three@0.170.0`, `webm-muxer@5.1.4`,
  `mp4-muxer@5.2.2`. `.mjs` files are renamed to `.js`. Each licence is
  recorded in `web/vendor/LICENSES.txt` and in `NOTICE`.
- Default visual style is `dark`. The look is an original take on the genre and
  must not copy any specific existing animation.
- Export sizes are exactly `1280x720` and `1920x1080`. Frame rate is 30 fps.
  Maximum script duration is `300` s.
- Protein actors fetch structure only through `GET /api/analyze?pdb=<ID>&chain=<C>`
  (same origin). One request is made per chain, and `pdb_data` is read from the
  response.
- World units: 1 unit = 1 nm for real structures (PDB Å ÷ 10). Cell-scale
  actors use stylised sizes, and the scene is labelled "not to scale" via
  `note`.
- **Bare-surface rule** (from `tests/test_ui_bare.py`): the studio editor opens
  empty. Examples load only when picked from the Examples menu. No
  pre-filled PDB IDs appear in `studio.html`.
- Every new `.js` file must pass `node --check`. Pure modules must import
  nothing from `render/`, `vendor/` or the DOM.
- Errors and warnings are `{line, message}`, with 1-based line numbers and
  plain-English messages.
- Commit messages end with the session's Co-Authored-By and Claude-Session
  lines.

## Review Focus

1. **Script edited mid-typing (half-finished lines).** The preview keeps the
   last good scene and shows errors. It never throws or blanks the canvas.
   Tested in Task 9 (`test_studio.py::test_broken_script_keeps_last_good_preview`).
2. **A PDB fetch that fails or names a missing chain** (offline, typo,
   `too_large`). That protein renders as a generic blob, and a warning names
   the line. Playback and export still work. Tested in Task 5
   (`structures.test.js` "falls back on 404 / missing chain").
3. **A reference to an actor declared later, or never declared**
   (`to=dnna@0.4`). Compile error naming the line, with a "did you mean"
   suggestion. Tested in Task 4 (`timeline.test.js` "unknown ref").
4. **Export in a browser without WebCodecs or the requested codec.** The
   Export menu lists only supported formats. If none are supported it is
   disabled with an explanation. Tested in Task 10 (`record.test.js`
   "supportedFormats with no VideoEncoder").
5. **Seeking backwards or jumping to any time.** `stateAt` is pure, so
   `stateAt(t)` is identical whether reached by play, seek or export. Tested in
   Task 4 (`timeline.test.js` "stateAt is order-independent").

---

## File map

```
web/studio.html                      page shell (editor | preview | timeline)
web/studio.css                       studio chrome; reuses web/style.css tokens
web/studio/package.json              {"type":"module","private":true} so Node treats .js as ESM
web/studio/schema.js                 actor types, verbs, header keys, param specs, suggest()
web/studio/parse.js                  text -> Scene (syntax + schema validation)
web/studio/anchors.js                vector/euler helpers, actorExtent(), anchorPoint()
web/studio/timeline.js               compile(Scene) -> Timeline; stateAt(); ease()
web/studio/pdbcoarse.js              caCoords(), coarseBeads()
web/studio/structures.js             loadStructures(scene, fetchFn) -> Map(actorName -> Beads|null)
web/studio/render/materials.js       glowMaterial()
web/studio/render/molecular.js       helix, strand, protein, point, arrow views
web/studio/render/cellular.js        ring, chromosome, fibre, cell views
web/studio/render/stage.js           Stage: scene graph, camera, apply(FrameState), project()
web/studio/render/overlay.js         drawOverlay(): captions, labels, note, credits
web/studio/record.js                 supportedFormats(), exportVideo()
web/studio/app.js                    page wiring
web/studio/examples/transcription-basic.txt
web/studio/examples/nondisjunction-meiosis-1.txt
web/studio/test/*.test.js            node --test suites
web/vendor/{three.module.min.js, MarchingCubes.js, webm-muxer.js, mp4-muxer.js, LICENSES.txt}
server.py                            + ".txt" content type
tests/test_studio_static.py          static serving + nav link + bare editor
e2e/test_studio.py                   Playwright smoke
.github/workflows/ci.yml             + node --check web/studio/**, node --test
```

## Shared data shapes (the contract between tasks)

```js
// parse.js
ParseResult = { scene: Scene|null, errors: Diag[], warnings: Diag[] }
Diag   = { line: number, message: string }
Scene  = { title: string, size: [w,h], duration: number|null, style: 'dark'|'light',
           note: string|null, credits: boolean,
           actors: ActorDecl[], actions: ActionDecl[] }
ActorDecl  = { name, type, params: {[k]: Value}, line }
ActionDecl = { t0, t1, target, verb, args: {[k]: Value}, text: string|null, line }
                       // target: actor name | 'camera' | 'caption'; text: quoted caption text
Value = number | string | boolean | '#rrggbb' | Vec3 | Ref | {dir: Dir}
Vec3 = [x,y,z]
Ref  = { ref: string, u?: number, part?: string }   // "dna@0.45" | "chr21a.kinetochore" | "sep"
Dir  = 'left'|'right'|'top'|'bottom'|'front'|'back'

// timeline.js
Timeline   = { duration, size, style, note, creditsText: string|null,
               actors: {[name]: {type, params, line}}, tracks: Track[],
               captions: {t0,t1,text}[], cameraTracks: Track[], warnings: Diag[] }
Track      = { actor, prop, t0, t1, from, to, ease, line, dynamic?: Ref }
FrameState = { t, actors: {[name]: ActorState}, camera: {pos: Vec3, target: Vec3, fov: number},
               captions: string[], creditsVisible: boolean }
ActorState = { visible: boolean, opacity: number, pos: Vec3, rot: Vec3 /*deg XYZ*/,
               scale: number, color: '#rrggbb', glow: number,
               p: {[k]: number} /* type-specific, see schema STATE_PROPS */ }

// pdbcoarse.js / structures.js
Beads = { beads: [x,y,z,r][], radius: number }   // centred on origin, nm units
```

---

### Task 1: Scaffold, vendored libraries, serving, CI

**Files:**
- Create: `web/vendor/three.module.min.js`, `web/vendor/MarchingCubes.js`, `web/vendor/webm-muxer.js`, `web/vendor/mp4-muxer.js`, `web/vendor/LICENSES.txt`
- Create: `web/studio.html`, `web/studio.css`, `web/studio/package.json`, `web/studio/app.js` (stub)
- Modify: `server.py:104-111` (`_CONTENT_TYPES`), `NOTICE`, `.github/workflows/ci.yml`
- Modify: the sitenav in `web/index.html`, `web/help.html`, `web/api.html`, `web/privacy.html` and `web/terms.html`. Add `<a href="/studio.html">Studio</a>` after Workbench.
- Test: `tests/test_studio_static.py`

**Interfaces:**
- Produces: `/studio.html` loads `studio/app.js` as `<script type="module">`. It has these element IDs: `#editor` (textarea), `#diagnostics`, `#preview` (canvas container), `#overlay` (2D canvas stacked on the WebGL canvas), `#play`, `#scrub` (range input), `#clock`, `#examples` (select), `#open` (file input), `#save`, `#export` (select or button group), `#exportStatus`.

- [ ] **Step 1: Write failing tests** in `tests/test_studio_static.py`, following the server-boot pattern in `tests/test_server_integration.py`.
  - `test_studio_page_served`: GET `/studio.html` returns 200 with `text/html`, and the body contains `type="module"` and `studio/app.js`.
  - `test_example_txt_served_as_text`: GET `/studio/examples/transcription-basic.txt` returns 200 with `Content-Type` starting `text/plain`. Write a placeholder example file in this step; Task 11 replaces it.
  - `test_vendor_three_has_no_bare_imports`: in each `web/vendor/*.js`, no `from 'three'` or `from "three"`.
  - `test_studio_linked_from_every_page`: each of the five HTML pages contains `href="/studio.html"`.
  - `test_editor_starts_empty`: in `studio.html`, the `<textarea id="editor">` element has no text content, and the file contains no 4-character PDB-ID-looking token from `tests/test_ui_bare.py`'s demo list.
- [ ] **Step 2: Run** `python -m unittest tests.test_studio_static -v`. Expected: FAIL (404s, missing files).
- [ ] **Step 3: Vendor the libraries.**
  - Download them with `npm pack`.
  - Copy `build/three.module.min.js` and `examples/jsm/objects/MarchingCubes.js`. Rewrite MarchingCubes' import source to `'./three.module.min.js'`.
  - Copy `build/webm-muxer.mjs` to `webm-muxer.js`, and `build/mp4-muxer.mjs` to `mp4-muxer.js`.
  - In `LICENSES.txt`, add the three MIT licence texts with versions.
  - In `NOTICE` under SOFTWARE, add one entry per library, noting "vendored in web/vendor/".
- [ ] **Step 4:** Add `".txt": "text/plain; charset=utf-8"` to `_CONTENT_TYPES`. Write the `studio.html` and `studio.css` shell with the IDs above, and the sitenav with `Studio` marked `class="here"`. Make `app.js` a stub that imports `../vendor/three.module.min.js` and logs nothing. Add the nav links to the other pages.
- [ ] **Step 5:** Add two steps to `ci.yml` after "Syntax-check the front end":
  - `find web/studio -name '*.js' -print0 | xargs -0 -n1 node --check`
  - `node --test web/studio/test/`

  Add an empty `web/studio/test/smoke.test.js` with one passing test, so the runner has a target.
- [ ] **Step 6: Run** `python -m unittest discover -s tests` and `node --test web/studio/test/`. Expected: all PASS, including the existing `test_ui_bare` and `test_security_headers_on_static`.
- [ ] **Step 7: Commit** with the message "Scaffold Studio page, vendor Three.js and muxers".

---

### Task 2: Script grammar and parser (syntax)

**Files:**
- Create: `web/studio/parse.js`
- Test: `web/studio/test/parse.test.js`

**Interfaces:**
- Produces: `parse(text: string) → ParseResult`. In this task, values are typed by syntax only. Schema checks come in Task 3, which calls `validate(scene)` from inside `parse`.

Grammar (line-based; `#` starts a comment, except inside quotes; blank lines are ignored):

```
header   := 'title' STRING | 'size' INT 'x' INT | 'duration' TIME | 'style' WORD
          | 'note' STRING | 'credits' ('on'|'off')
actor    := 'actor' NAME '=' TYPE (KEY '=' VALUE)*
action   := TIMESPEC TARGET VERB [OBJECT] (KEY '=' VALUE)*
           -- OBJECT is one bare ref after the verb, stored as args.object
           --   ("fibre1 pull chr21a to=poleL", "camera zoom ring"); Task 3 maps it to the verb's
           --   object arg: pull->target, zoom->to, follow->ref, orbit->around, pair->with,
           --   attach->to, camera cut->to. An OBJECT on any other verb is an error.
           -- caption form: TIMESPEC 'caption' STRING (KEY=VALUE)*
TIMESPEC := TIME | TIME '-' TIME          TIME := NUMBER ['s']      (e.g. 6, 6s, 0-4s, 0s-4s, 2.5-3.25s)
VALUE    := NUMBER | STRING | '#'HEX6 | COLORNAME | X','Y','Z | NAME'@'NUMBER | NAME'.'PART | NAME | DIR | yes|no|on|off
NAME     := [A-Za-z][A-Za-z0-9_-]*
```

Rules:
- A single-time action has `t1 = t0 + defaultDuration(verb)`. The table lives in `schema.js`, and parse gets it from `validate` (Task 3). In this task, store `t1 = null`.
- `TIME '-' TIME` with `t1 < t0` is an error: `"end time is before start time"`.
- Colour names map to hex: pink `#ff5fa2`, red `#ff4d4d`, orange `#ff9a3c`, gold `#ffc94d`, yellow `#fff06b`, green `#5fe07a`, teal `#3fd6c6`, blue `#4da3ff`, violet `#a77bff`, white `#f2f2f2`, grey `#9aa0a6`.

- [ ] **Step 1: Write the failing tests** (`node:test`, `node:assert/strict`).
  - "parses header": `title "X"\nsize 1280x720\nduration 45s\nstyle dark\nnote "n"\ncredits off` gives `scene.title === 'X'`, `size` deep-equals `[1280,720]`, `duration === 45`, `credits === false`.
  - "parses actor with typed values": `actor sep = protein from=7NJ1 chain=A color=gold at=1,2,3` gives params `{from:'7NJ1', chain:'A', color:'#ffc94d', at:[1,2,3]}`, `line === 1`.
  - "parses refs": `to=dna@0.45` gives `{ref:'dna', u:0.45}`. `to=chr21a.kinetochore` gives `{ref:'chr21a', part:'kinetochore'}`. `from=left` gives `{dir:'left'}`.
  - "parses ranges and instants": `0-4s sep enter from=left` gives `t0 0, t1 4, verb 'enter'`. `6s ring open` gives `t0 6, t1 null`.
  - "parses positional object": `14-20s fibre1 pull chr21a to=poleL` gives `verb 'pull'`, `args.object` deep-equal to `{ref:'chr21a'}`, and `args.to` deep-equal to `{ref:'poleL'}`.
  - "caption with text": `6s caption "Securin released" for=3s` gives `target 'caption'`, `text 'Securin released'`, `args.for === 3`.
  - "comments and quotes": `title "A # B" # c` gives title `'A # B'`.
  - "reversed range errors with line": line 3 `5-2s x move` gives `errors[0]` deep-equals `{line:3, message:'end time is before start time'}`.
  - "unterminated string": the error message contains `missing closing "`.
  - "never throws": for each of 200 random truncations of the nondisjunction sketch in the spec, `parse` returns an object and does not throw.
- [ ] **Step 2: Run** `node --test web/studio/test/parse.test.js`. Expected: FAIL (module not found).
- [ ] **Step 3: Implement `parse`** as a per-line tokenizer (quoted strings first, then whitespace split), then line classification by first token. Collect errors instead of throwing. `scene` is `null` only when there are errors.
- [ ] **Step 4: Run tests.** Expected: PASS.
- [ ] **Step 5: Commit** with the message "Studio: scene-script parser".

---

### Task 3: Schema validation and "did you mean"

**Files:**
- Create: `web/studio/schema.js`
- Modify: `web/studio/parse.js` (call `validate` before returning)
- Test: `web/studio/test/schema.test.js`

**Interfaces:**
- Produces:
  - `ACTOR_TYPES: {[type]: {params: {[key]: ParamSpec}, parts: string[], stateProps: {[prop]: number /*default*/}}}`
  - `VERBS: {[verb]: {appliesTo: string[] /*types or '*'*/, args: {[key]: ParamSpec}, defaultDuration: number}}`
  - `CAMERA_VERBS`
  - `suggest(word: string, candidates: string[]) → string|null`, using Levenshtein distance ≤ 2 and the closest match.
  - `validate(scene: Scene) → {errors: Diag[], warnings: Diag[]}`. It fills `t1` for instants and coerces defaults.
- `ParamSpec = {kind: 'number'|'string'|'color'|'vec3'|'ref'|'dir'|'refOrVec'|'refOrDir'|'bool'|'time', default?, min?, max?, required?}`

Actor types and type-specific state props (defaults in parentheses):

| type | params (beyond common `at`, `color`, `scale`, `rot`, `label`, `hidden`) | parts | stateProps |
|---|---|---|---|
| helix | `turns` (20, 1–200), `curve` (0, 0–1) | `start`, `end` | `unwind0` (0), `unwind1` (0), `cut` (0), `cutAt` (0.5) |
| strand | `length` (10), `along` ref | `tip` | `grow` (1) |
| protein | `from` (PDB ID, optional), `chain` (string; comma list allowed), `size` (6) | — | — |
| ring | `around` ref, `radius` (3) | — | `open` (0), `cleave` (0) |
| chromosome | `sisters` (2, 1–2), `length` (12) | `centromere`, `kinetochore`, `top`, `bottom` | `separate` (0), `condense` (1) |
| fibre | `from` refOrVec, `to` refOrVec | `tip` | `grow` (1) |
| cell | `radius` (40), `nucleus` (bool, yes) | — | `pinch` (0), `split` (0), `envelope` (1) |
| point | `show` (bool, no) | — | — |
| label | `text` (required), `on` ref | — | — |
| arrow | `from` ref, `to` ref | — | — |

Verbs (`appliesTo`, args, defaultDuration in s):
- On any actor: `enter` (from refOrDir, to refOrVec; 1.5), `leave` (to refOrDir; 1.5), `move` (to refOrVec; 1.5), `attach` (to ref; 0.5), `detach` (0), `rotate` (by number deg, axis x|y|z; 1.5), `scale` (to number; 1), `fade` (to number 0–1; 1), `color` (to color; 1), `glow` (to number 0–1; 1), `show` (0.5), `hide` (0.5).
- Helix: `unwind` (from 0–1, to 0–1; 2), `rewind` (2), `cut` (at 0–1, gap number default 4; 1.5).
- Strand: `grow` (to 0–1; 2), `peel` (1.5).
- Ring: `open`, `close`, `cleave` (each 1).
- Chromosome: `separate` (gap default 6; 2), `condense` (to 0–1; 2), `pair` (with ref; 2).
- Fibre: `grow` (2), `shrink` (2), `pull` (target ref, to refOrVec; 3).
- Cell: `pinch` (to 0–1; 2), `split` (gap default 50; 3).
- Camera: `wide` (1.5), `zoom` (to ref, distance number; 2), `orbit` (around ref, degrees default 360; duration from range, else 6), `follow` (ref; duration from range, else 4), `cut` (to ref; 0).
- Caption: arg `for` (time, default 3). If a range is given, the range wins.

Validation rules:
- Unknown header key, actor type, verb or param key → error ending in `— did you mean "x"?` when `suggest` returns a match.
- A verb that doesn't apply to the target's type → error: `"<verb>" does not apply to <type> "<name>"`.
- A target name that isn't declared → error, with a suggestion from declared names.
- Every `Ref` in actor params or action args must name a declared actor
  (declaration order does not matter). An undeclared name → error with a
  suggestion from the declared names. A `.part` must be in that type's `parts`,
  and a `@u` is only allowed on `helix`, with 0 ≤ u ≤ 1.
- Duplicate actor name → error naming the earlier line.
- Number out of range → error with the allowed range.
- `duration > 300` → error `"maximum duration is 300s"`. `size` not in {1280x720, 1920x1080} → error listing both. Style not in dark/light → error.
- A protein `from` must match `/^[0-9][A-Za-z0-9]{3}$/`.

- [ ] **Step 1: Write failing tests.**
  - "suggests close words": `suggest('pul', ['pull','push'])` returns `'pull'`, and `suggest('zzz', ['pull'])` returns `null`.
  - "unknown verb": `parse('actor f = fibre\n1s f pul')` gives an error with `{line:2}` whose message ends with `did you mean "pull"?`.
  - "verb not for type": `actor d = helix\n1s d open` gives the message `"open" does not apply to helix "d"`.
  - "undeclared target with suggestion": `actor dna = helix\n0-2s dnna unwind from=0.2 to=0.4` gives an error containing `did you mean "dna"?`.
  - "instant gets default duration": `actor r = ring\n6s r open` gives `actions[0].t1 === 7`.
  - "duration cap": `duration 301s` gives an error containing `300s`.
  - "bad size": `size 800x600` gives an error mentioning `1280x720` and `1920x1080`.
  - "duplicate actor": the error message contains `already declared on line 1`.
  - "every verb in VERBS has defaultDuration and appliesTo": an iteration assertion.
- [ ] **Step 2: Run** `node --test web/studio/test/schema.test.js`. Expected: FAIL.
- [ ] **Step 3: Implement `schema.js`** with the tables above, and call `validate` from `parse`.
- [ ] **Step 4: Run** `node --test web/studio/test/`. Expected: PASS, and parse tests still pass.
- [ ] **Step 5: Commit** with the message "Studio: schema validation with suggestions".

---

### Task 4: Timeline compiler and `stateAt`

**Files:**
- Create: `web/studio/anchors.js`, `web/studio/timeline.js`
- Test: `web/studio/test/timeline.test.js`, `web/studio/test/anchors.test.js`

**Interfaces:**
- Consumes: `Scene` and `ACTOR_TYPES`/`VERBS` from Tasks 2–3.
- Produces:
  - `anchors.js`:
    - `rotate(v: Vec3, rotDeg: Vec3) → Vec3`, applying XYZ Euler order.
    - `actorExtent(type, params) → number`, the bounding radius in world units.
    - `anchorPoint(decl: ActorDecl, state: ActorState, ref: Ref) → Vec3`.
  - `timeline.js`:
    - `ease(name: 'smooth'|'linear'|'in'|'out', x: number) → number`. `smooth` is cubic in-out.
    - `compile(scene: Scene) → {timeline: Timeline|null, errors, warnings}`.
    - `stateAt(timeline: Timeline, t: number) → FrameState`.

Geometry conventions (anchors and renderer must agree):
- **Helix:** the axis lies along local +x, centred on `pos`. Length is `turns × 3.4` (rise 0.34 nm/bp × 10 bp per turn), and radius is 1. `dna@u` is `pos + rotate([(u − 0.5) × L, 0, 0]) × scale`.
- **Chromosome:** the long axis is local y with length `length`. `centromere` is `pos`. `kinetochore` is `pos + rotate([−1.2, 0, 0])`, the face pointing to local −x. `top` and `bottom` are ±length/2 on y.
- **Fibre:** `from` and `to` resolve dynamically every frame. `tip` is `lerp(from, to, grow)`.
- **Strand:** if `along` is set, it lies along that helix. `tip` is the helix point at `u0 + grow × (u1 − u0)` of the target helix's current unwind range.
- `actorExtent`:
  - helix: `L/2`
  - chromosome: `length/2`
  - cell: `radius`
  - ring: `radius`
  - protein: `size/2`, until Beads radius is known, which the stage may override
  - point, label and arrow: `1`
- `Dir` offsets for `enter`/`leave` are 60 units along: left −x, right +x, top +y, bottom −y, front +z, back −z.

Compile semantics:
- Process actions sorted by `(t0, line)`. Each verb expands to one or more `Track`s on props:
  - `pos`, `rot`, `scale`, `opacity`, `color`, `glow`, `visible`
  - `p.<stateProp>`, for type-specific state
- A track's `from` is the prop's value at `t0` given the tracks compiled so far. `to` comes from the verb.
- Refs in `to=` or `from=` resolve to points **at compile time**, at the track's t1 (`to`) or t0 (`from`), using `stateAt` over the tracks compiled so far.
- `attach` and `follow`, and fibre endpoints, are **dynamic**. They create a track with `dynamic: Ref`, which `stateAt` resolves every frame. Each frame resolves actors in dependency order. A cycle is a compile error: `"attach cycle: a → b → a"`.
- Visibility:
  - If an actor's first action is `enter` or `show`, it is invisible before that action. Otherwise it is visible from t = 0, unless `hidden=yes`.
  - `leave` ends with `visible: false` and `opacity: 0`.
- Overlap: two tracks on the same `(actor, prop)` whose intervals overlap → the later `(t0, line)` wins from its t0. A warning is added: `line X overrides line Y (<prop> of "<actor>")`.
- Duration:
  - If `scene.duration` is set and an action ends after it → error `action ends at <t1>s but duration is <d>s`.
  - If `scene.duration` is null → `duration = max(t1) + 1`, plus 4 more when credits are on.
- Credits text:
  - If `credits` is on and at least one protein has `from`: `"Structures: " + ", ".join("<ID> chain <C>" + (label ? " (<label>)" : ""))` in declaration order.
  - Append `" · Rendered with SnaCleX Studio"`.
  - It is visible during `[duration − 4, duration]`.
- Camera:
  - With no camera actions, the camera is `wide` over all actors visible at t = 0: target at the centroid of their positions, distance `1.6 × max(extent + |pos − centroid|)`, viewing from +z with slight elevation `[0, 0.25, 1]` normalised. `fov` is 35.
  - `zoom to=ref distance=d`: the default `d` is `3 × extent` of the ref's actor.
  - `orbit` rotates the camera position about the target's y axis by `degrees` over the interval.
  - `follow` makes the target dynamic.
  - `cut` is instant.

- [ ] **Step 1: Write failing tests.**
  - anchors: "helix anchor": a helix with turns 10 at `[0,0,0]` with no rotation gives `dna@0` = `[-17,0,0]` and `dna@1` = `[17,0,0]` (within 1e-9). With `rot:[0,90,0]`, `dna@1` ≈ `[0,0,-17]`.
  - "ease endpoints": for each easing, `ease(e,0) === 0` and `ease(e,1) === 1`, and `ease('smooth',0.5) === 0.5`.
  - "move interpolates position": `actor p = protein at=0,0,0\n0-2s p move to=10,0,0` gives `stateAt(tl,1).actors.p.pos[0]` = 5 (smooth at the midpoint). At `t=3` it is 10.
  - "enter hides before start": `actor p = protein\n2-4s p enter from=left to=0,0,0` gives `visible === false` at `t=1`. At `t=4`, `pos` deep-equals `[0,0,0]` and `opacity === 1`.
  - "ref resolves to helix point": `actor dna = helix turns=10\nactor p = protein\n0-1s p move to=dna@1` gives `p.pos` ≈ `[17,0,0]` at `t=1`.
  - "attach follows a moving target": p is attached at 0–0.5s to q, then `1-2s q move to=0,10,0`. At `t=2`, `p.pos[1]` is 10 plus the original offset.
  - "unknown ref" (Review Focus 3): `0-1s p move to=dnna@0.4` gives an error containing `did you mean "dna"?`. (Parse/validate may catch this first. Either way, assert the combined `errors`.)
  - "attach cycle": the error message starts `attach cycle`.
  - "overlap warning": two overlapping moves on p → one warning matching `/overrides line \d+ \(pos of "p"\)/`.
  - "duration overflow": `duration 5s` with an action ending at 6 → the error message is `action ends at 6s but duration is 5s`.
  - "credits text": `actor sep = protein from=7NJ1 chain=A label="separase"` gives `creditsText` that starts `Structures: 7NJ1 chain A (separase)`, and `stateAt(tl, tl.duration - 1).creditsVisible === true`.
  - "stateAt is order-independent" (Review Focus 5): for a 6-action scene, `JSON.stringify(stateAt(tl,3))` is the same before and after calling `stateAt` at t = 9, 0 and 5.
  - "caption visibility": `2s caption "Hi" for=3s` gives `captions` deep-equal to `['Hi']` at t = 4, and `[]` at t = 5.5.
- [ ] **Step 2: Run** `node --test web/studio/test/`. Expected: FAIL.
- [ ] **Step 3: Implement `anchors.js`.**
- [ ] **Step 4: Implement `timeline.js`.** `stateAt` must not mutate `timeline`. Return fresh objects.
- [ ] **Step 5: Run tests.** Expected: PASS.
- [ ] **Step 6: Commit** with the message "Studio: timeline compiler and pure stateAt".

---

### Task 5: Real-structure proteins (coarse beads)

**Files:**
- Create: `web/studio/pdbcoarse.js`, `web/studio/structures.js`
- Test: `web/studio/test/pdbcoarse.test.js`, `web/studio/test/structures.test.js`

**Interfaces:**
- Produces:
  - `caCoords(pdbText: string, chain: string) → Vec3[]`. It reads `ATOM` records whose atom name is `CA` in the given chain, using PDB fixed columns, and divides by 10 to get nm.
  - `coarseBeads(coords: Vec3[], perBead = 4) → Beads`. It averages each run of `perBead` consecutive residues into one bead of radius 0.9, centres everything on the centroid, and sets `radius` to the max distance plus 0.9.
  - `loadStructures(scene: Scene, fetchFn = fetch) → Promise<{beads: Map<string, Beads|null>, warnings: Diag[]}>`. For each protein actor with `from`, and each chain in its comma list, it calls `GET /api/analyze?pdb=<ID>&chain=<C>`. It concatenates the CA coordinates across chains before `coarseBeads`. A module-level cache is keyed by `ID:chains`.
  - On a non-OK response, JSON with `error` or `too_large`, zero CA atoms, or a thrown fetch: the map value is `null`, plus the warning `structure <ID> chain <C> could not be loaded (<reason>); showing a generic shape`. The line is the actor's line.

- [ ] **Step 1: Write failing tests.**
  - "caCoords reads CA only": a 3-residue fixture string with N, CA and C atoms in chain A, plus one CA in chain B, gives 3 coordinates for chain A, the first equal to the fixture's CA ÷ 10.
  - "coarseBeads centres and sizes": 8 collinear coordinates and `perBead` 4 give 2 beads with a centroid of about `[0,0,0]` and `radius > 0.9`.
  - "loads via /api/analyze with chain": a fake `fetchFn` records the URL and returns `{ok:true, json:()=>({pdb_data: FIXTURE})}`. The URL is `/api/analyze?pdb=7NJ1&chain=A`, and the map has a `Beads` value.
  - "falls back on 404 / missing chain" (Review Focus 2): with a fake returning `{ok:false, status:404}`, the value is `null` and there is one warning on the actor's line containing `showing a generic shape`. The same holds for JSON `{too_large:true}`.
  - "multi-chain concatenates": `chain=A,B` makes two fetches, and the bead count reflects both chains.
- [ ] **Step 2: Run tests.** Expected: FAIL.
- [ ] **Step 3: Implement both modules.**
- [ ] **Step 4: Run tests.** Expected: PASS.
- [ ] **Step 5: Commit** with the message "Studio: coarse protein beads from SnaCleX structures".

---

### Task 6: Stage, materials and molecular actor views

**Files:**
- Create: `web/studio/render/materials.js`, `web/studio/render/molecular.js`, `web/studio/render/stage.js`
- Modify: `web/studio/app.js`. Add a temporary dev hook that renders a hard-coded scene so the stage can be seen; Task 9 replaces it.

**Interfaces:**
- Consumes: `Timeline`, `FrameState` and `anchors.js` geometry conventions (Task 4); `Beads` (Task 5).
- Produces:
  - `glowMaterial(color: string, style: 'dark'|'light') → THREE.Material`. Use `MeshStandardMaterial` with a Fresnel rim emissive term added via `onBeforeCompile`. Opacity and colour can be changed per frame.
  - `ActorView = {object: THREE.Object3D, update(s: ActorState, f: FrameState, tl: Timeline): void, dispose(): void}`.
  - `makeHelix(decl)`, `makeStrand(decl)`, `makeProtein(decl, beads: Beads|null)`, `makePoint(decl)`, `makeArrow(decl)`. Each returns an `ActorView`.
  - `class Stage`:
    - `constructor(canvas: HTMLCanvasElement, {style})`
    - `build(tl: Timeline, beads: Map): void`
    - `apply(f: FrameState): void`
    - `render(): void`
    - `setSize(w, h): void`
    - `project(v: Vec3) → [x, y] | null`, giving canvas pixels, or null when behind the camera
    - `dispose(): void`
  - The stage registers its builders in a `VIEW_BUILDERS` map keyed by actor type. Task 7 adds the cellular ones to the same map.
  - The canvas is created with `preserveDrawingBuffer: true`, so the export can read frames.

Visual spec:
- **Dark style:** background `#05060a`, exponential fog `0.006`, one hemisphere light, and one key directional light.
- **Light style:** background `#f4f1ea` and no fog.
- **Helix:** two backbone tubes (`TubeGeometry` along helical curves, radius 0.18) and base-pair rods every 0.34 units.
  - `unwind0..unwind1` opens a bubble in that `u` range: the strands separate to radius 2.2 with smoothstep falloff at the edges, and base pairs inside the bubble are hidden.
  - `cut` > 0 splits the helix at `cutAt` and translates the halves apart by `gap × cut`. The `gap` value is stored on the track.
  - Geometry is rebuilt only when `unwind0`, `unwind1` or `cut` changes by more than 1e-3.
- **Protein:**
  - With beads: a `MarchingCubes` field (resolution 48) with one ball per bead, scaled to `Beads.radius`. This gives a smooth blob.
  - Without beads: a generic lumpy blob from 5 seeded balls in a `size`-diameter cluster, seeded by an FNV-1a hash of the actor name, so it is stable across reloads.
- **Strand:** a single tube along the bubble of its `along` helix, drawn from `u0` to the grow point, radius 0.15.
- **Point:** invisible unless `show=yes`, in which case it is a small glowing sphere. **Arrow:** `ArrowHelper`-like cone and shaft between its two refs.

- [ ] **Step 1:** With the temporary dev hook, render this scene: a helix (turns 20), a protein (no `from`) moving onto `dna@0.5` over 0–3s, and the helix unwinding 0.45–0.55 over 3–5s. Do it through `compile` → `Stage.build` → per-frame `apply(stateAt(...))`.
- [ ] **Step 2: Verify visually.**
  - Serve the app with `python server.py`.
  - Use a Playwright script in the scratchpad to load `/studio.html?dev=1`, seek to t = 0, 2.5 and 5, and screenshot each (deviceScaleFactor 1).
  - Check: dark background, a pink double helix, the blob reaching the helix by t = 2.5, and a visible bubble at t = 5.
  - Check there are no console errors.
- [ ] **Step 3:** Run `node --check` on the new files. Expected: OK.
- [ ] **Step 4: Commit** with the message "Studio: stage, glow material, molecular actor views".

---

### Task 7: Cellular actor views

**Files:**
- Create: `web/studio/render/cellular.js`
- Modify: `web/studio/render/stage.js`. Register `ring`, `chromosome`, `fibre` and `cell` in `VIEW_BUILDERS`.

**Interfaces:**
- Consumes: `ActorView`, `glowMaterial` and `VIEW_BUILDERS` (Task 6); the anchor conventions (Task 4).
- Produces: `makeRing(decl)`, `makeChromosome(decl)`, `makeFibre(decl)`, `makeCell(decl)`, each returning an `ActorView`.

Visual spec:
- **Ring:** a torus of the given `radius` (tube 0.35), placed perpendicular to the `around` target's axis.
  - `open` 0–1 opens a gap that grows to 90° (the torus arc shrinks from 2π to 1.5π).
  - `cleave` 0–1 splits it into two arcs that drift apart by 3 units and fade by 50%.
- **Chromosome:** each sister chromatid is a capsule of length `length` with a pinched waist at the centromere. Draw it with `LatheGeometry`, with the waist at 55% radius.
  - Two sisters sit side by side, ±1.1 on local x, joined at the centromere.
  - `separate` moves the sisters apart by `gap × separate` on local x.
  - `condense` scales radius from 0.6 to 1.0 and length from 1.4× to 1.0×.
  - A small glowing disc marks the `kinetochore` part.
- **Fibre:** a thin tube (radius 0.12) from the resolved `from` to `lerp(from, to, grow)`. Colour defaults to `#9fe3ff` when unset.
- **Cell:**
  - A translucent sphere (opacity 0.12, glow rim) of `radius`.
  - The nucleus is an inner sphere at 0.45 × radius whose opacity is `envelope`.
  - `pinch` squeezes the equator inward toward 0.2 × radius using a vertex-shader or morph displacement.
  - `split` replaces the cell with two hemispherical cells that move apart by `gap × split`.

- [ ] **Step 1:** Extend the dev hook. Show a chromosome pair, with a cohesin ring around them that opens over 1–2s; a fibre from a `point` at `-30,0,0` growing to the kinetochore; and a cell pinching over 3–5s.
- [ ] **Step 2: Verify visually.** Using the same Playwright screenshot routine at t = 0, 1.5, 3 and 5, confirm each effect is visible and there are no console errors.
- [ ] **Step 3: Commit** with the message "Studio: ring, chromosome, fibre and cell views".

---

### Task 8: Camera, overlay (captions, labels, note, credits)

**Files:**
- Create: `web/studio/render/overlay.js`
- Modify: `web/studio/render/stage.js`. `apply` sets the camera from `f.camera`.
- Test: `web/studio/test/overlay.test.js` (pure layout helper only)

**Interfaces:**
- Produces:
  - `layoutCaption(text: string, width: number, measure: (s) => number, maxWidthFrac = 0.8) → string[]`. It word-wraps to lines that fit within `width × maxWidthFrac`. It is pure, so it can be tested.
  - `drawOverlay(ctx: CanvasRenderingContext2D, f: FrameState, tl: Timeline, stage: Stage, {width, height, style}) → void`. It draws:
    - **Captions:** bottom-centred, font size `height × 0.038`, white on a semi-transparent black band in dark style (inverted in light style).
    - **Labels:** for every `label` actor and every actor with a `label` param, when that actor is visible. The text sits at `stage.project(anchor)`, offset up by 1.5 × extent, with a thin leader line.
    - **Note:** top-left, small, 60% opacity, for the whole video.
    - **Credits:** a centred card when `creditsVisible`, with the title on the first line and `creditsText` below.

- [ ] **Step 1: Write failing tests.**
  - "wraps long captions": with `measure = s => s.length * 10` and width 300, `layoutCaption('the quick brown fox jumps over', 300, measure)` gives lines that are each ≤ 240 px, and joining them with a space returns the original text.
  - "single short line unchanged": `layoutCaption('Hi', 1000, measure)` deep-equals `['Hi']`.
- [ ] **Step 2: Run** `node --test web/studio/test/overlay.test.js`. Expected: FAIL.
- [ ] **Step 3: Implement** `overlay.js` and the camera in `Stage.apply`. Use a `PerspectiveCamera` with `fov`, `position` and `lookAt(target)`.
- [ ] **Step 4:** Run the tests (expected: PASS). Take Playwright screenshots of a scene with a caption, a labelled protein and credits at the end, and confirm they read cleanly at 1280×720.
- [ ] **Step 5: Commit** with the message "Studio: camera moves and text overlay".

---

### Task 9: Studio page: editor, live preview, diagnostics, timeline, open/save/examples

**Files:**
- Modify: `web/studio/app.js` (remove the dev hook), `web/studio.html`, `web/studio.css`
- Test: `e2e/test_studio.py`. Follow `e2e/test_smoke.py`: boot the server, add a pre-seeded offline cache entry for a tiny 2-chain fixture structure `1ABC`, use Playwright, and skip if Playwright is absent.

**Interfaces:**
- Consumes: `parse`, `compile`, `stateAt`, `loadStructures`, `Stage` and `drawOverlay`.
- Produces: `window.__studio = {getTimeline(), seek(t), isPlaying()}` for E2E use only. It is attached only when `location.search` contains `test=1`.

Behaviour:
- On `input`, re-parse after 300 ms debounce.
  - If there are no errors: compile, load structures, and rebuild the stage. Keep the current time, clamped to the new duration.
  - If there are errors: **keep the last good stage and timeline**, and list the errors in `#diagnostics` as `Line N: message`. Clicking an entry moves the editor caret to that line.
  - Warnings are listed under the errors, in a muted style.
- Rendering: `requestAnimationFrame`. While playing, `t += dt` and it stops at the duration. `#scrub` sets `t` and pauses. `#clock` shows `m:ss.s / m:ss`.
- Preview size: the WebGL and overlay canvases fill `#preview` at 16:9, with device-pixel-ratio capped at 2.
- Examples: `#examples` lists "Transcription (basic)" and "Nondisjunction in Meiosis I", which load `studio/examples/*.txt` via fetch. Ask for confirmation if the editor has unsaved changes.
- Open: `<input type=file accept=".txt,text/plain">` puts the file text into the editor. Save: download the editor text as `<slugified title or "scene">.txt` via a Blob URL.
- The editor starts empty, with placeholder text describing the format only, not values. For example: "Write a scene script. See Help → Studio for the format."

- [ ] **Step 1: Write failing E2E tests.**
  - `test_example_loads_and_plays`: pick the transcription example. `#diagnostics` shows no `Line` errors. Click `#play`. After 1 s, `window.__studio.isPlaying()` is true. There are no `pageerror` events.
  - `test_broken_script_keeps_last_good_preview` (Review Focus 1): load an example, then append `\n3s ghost pul`. `#diagnostics` contains `did you mean`. `window.__studio.getTimeline()` still has the example's actors. There are no `pageerror` events.
  - `test_seek_updates_clock`: `window.__studio.seek(2)` makes `#clock` start with `0:02.0`.
  - `test_protein_from_fixture_structure`: the script `actor p = protein from=1ABC chain=A` produces no warning containing `generic shape`.
- [ ] **Step 2: Run** `python -m unittest e2e.test_studio -v`. Expected: FAIL.
- [ ] **Step 3: Implement** the behaviour above in `app.js`.
- [ ] **Step 4: Run** the E2E tests and the full unit suite. Expected: PASS.
- [ ] **Step 5: Commit** with the message "Studio: editor, live preview, diagnostics and timeline".

---

### Task 10: Video export

**Files:**
- Create: `web/studio/record.js`
- Modify: `web/studio/app.js` (Export menu and progress), `web/studio.html`
- Test: `web/studio/test/record.test.js`

**Interfaces:**
- Consumes: `Stage`, `drawOverlay`, `stateAt` and the vendored muxers.
- Produces:
  - `supportedFormats(env = globalThis) → Promise<('webm'|'mp4')[]>`.
    - `'webm'` is included if `env.VideoEncoder` exists and `isConfigSupported` accepts `vp09.00.40.08` or `vp8`.
    - `'mp4'` is included if it accepts `avc1.640028`.
    - If there is no `VideoEncoder`, it returns `[]`.
  - `exportVideo({timeline, stage, overlayCanvas, width, height, format, fps = 30, onProgress, signal}) → Promise<Blob>`.

Algorithm:
1. `stage.setSize(width, height)` and size the overlay canvas to match.
2. Create the muxer with an `ArrayBufferTarget`: `webm-muxer` for webm, or `mp4-muxer` with `fastStart: 'in-memory'` for mp4.
3. Configure `VideoEncoder` with `{codec, width, height, bitrate: width >= 1920 ? 12e6 : 6e6, framerate: fps}`.
4. For each `i` in `[0, ceil(duration × fps))`:
   - `f = stateAt(tl, i/fps)`, then `stage.apply(f)` and `stage.render()`, then `drawOverlay`.
   - Composite onto a 2D canvas: draw the WebGL canvas, then the overlay.
   - `new VideoFrame(composite, {timestamp: round(i × 1e6 / fps)})`, then `encode` with `keyFrame` every `2×fps` frames, then `close`.
   - Await `encoder.flush()` every 60 frames to bound memory.
   - Call `onProgress(i/total)`.
   - Check `signal.aborted` and throw `AbortError` if set.
5. Flush, finalize, and return `new Blob([target.buffer], {type: 'video/webm' | 'video/mp4'})`.
6. Restore the preview size in `finally`.

UI:
- The Export control offers `720p WebM`, `1080p WebM`, `720p MP4` and `1080p MP4`, filtered by `supportedFormats()`.
- If the list is empty, the control is disabled. `#exportStatus` reads "Video export needs a browser with WebCodecs (current Chrome, Edge, Firefox or Safari)."
- During export there is a progress bar and a Cancel button, and editing is disabled.
- The download is named `<slug>-<720p|1080p>.<ext>`.

- [ ] **Step 1: Write failing tests.**
  - "supportedFormats with no VideoEncoder" (Review Focus 4): `await supportedFormats({})` deep-equals `[]`.
  - "supportedFormats filters by codec": a fake `VideoEncoder.isConfigSupported` resolves `{supported: c.codec.startsWith('vp')}`, giving `['webm']`.
  - "frame timestamps are fixed-step": export a 1 s timeline at fps 30 with a fake stage and encoder that records `timestamp`s. There are 30 frames, the timestamps are `[0, 33333, 66667, …]`, and every 60th-frame `keyFrame` rule holds. Inject `VideoEncoder`, `VideoFrame` and the muxer via an optional `deps` param: `exportVideo(opts, deps = defaultDeps)`.
  - "abort stops export": aborting after 5 frames rejects with `name === 'AbortError'`.
- [ ] **Step 2: Run tests.** Expected: FAIL.
- [ ] **Step 3: Implement** `record.js` and the UI.
- [ ] **Step 4: Verify for real.**
  - In Playwright Chromium, load the transcription example and export 720p WebM via `page.expect_download()`.
  - Assert the file is larger than 50 kB.
  - If `ffprobe` exists, check its duration is within 0.2 s of the script duration. Otherwise check the `1A45DFA3` EBML magic bytes.
- [ ] **Step 5: Commit** with the message "Studio: deterministic WebM/MP4 export via WebCodecs".

---

### Task 11: The two example scripts

**Files:**
- Create: `web/studio/examples/transcription-basic.txt` and `web/studio/examples/nondisjunction-meiosis-1.txt`, replacing the Task 1 placeholder.
- Test: `web/studio/test/examples.test.js`

Content requirements:
- Both:
  - Have `title`, `style dark`, and `note "Stylised; not to scale"`.
  - Have `credits on`.
  - Are 30–90 s long.
  - Use captions for every step.
  - Use only 1–2 chains per protein, to keep fetches light.
- **Transcription (basic):**
  - Pol II is `from=1I6H chain=A,B` with `label="RNA polymerase II"`.
  - It binds the promoter, unwinds about 10% of the helix (the bubble), and moves along with a growing `strand` (mRNA, colour `#4da3ff`).
  - The bubble moves with it, using successive `unwind` actions.
  - Then the polymerase leaves and the helix rewinds.
- **Nondisjunction in Meiosis I** (the trisomy 21 story):
  1. Wide view of a `cell` with a nucleus. The envelope fades.
  2. The chromosome 21 homolog pair (two `chromosome` actors, maternal and paternal colours, `label`s) pairs at the centre. The camera zooms to the pair.
  3. A cohesin `ring` from `6WG3 chain=A,B` is shown on the chromatid arms. The ring is a stylised actor, and the real SMC1/SMC3 blob is shown alongside with `label="cohesin (SMC1/SMC3)"`.
  4. Separase (`7NJ1 chain=A`) is held by securin (`7NJ1 chain=B`). Securin leaves. Caption: "Securin is degraded; separase is released."
  5. Separase moves to the ring, and the arm cohesin is cleaved. Centromeric cohesin stays in Meiosis I, so a caption says so.
  6. Normal case, shown first: fibres from the two poles attach to each homolog's kinetochore, and the homologs pull to opposite poles. Then a `camera cut` to the error case.
  7. Error case: both homologs' kinetochores attach to fibres from the same pole. Both chr21 homologs are pulled to one side. Caption: "Nondisjunction: both homologs go to the same pole."
  8. The cell pinches and splits. Caption: "One gamete receives two copies (n+1), the other none (n−1)."
  9. Closing caption: "Fertilised by a normal gamete (n), the n+1 gamete gives a zygote with three copies of chromosome 21 (2n+1): trisomy 21."
- Scientific check: captions must match standard textbook accounts of Meiosis I nondisjunction (homologs fail to separate). They must not claim a single molecular cause is the established mechanism. Use "can", for example "Cohesin loss or mis-attachment can cause…".

- [ ] **Step 1: Write the failing test.** "examples parse and compile cleanly": for each file in `web/studio/examples/`, read it with `fs`, then `parse` gives `errors.length === 0`, `compile` gives `errors.length === 0` and `warnings.length === 0`, and `timeline.duration` is between 30 and 90.
- [ ] **Step 2: Run it.** Expected: FAIL (the placeholder is invalid or too short).
- [ ] **Step 3: Write both scripts.**
- [ ] **Step 4: Run** `node --test web/studio/test/`. Expected: PASS.
- [ ] **Step 5: Verify visually.** Export a 720p WebM of each example (Task 10 routine), extract frames at 25%, 50% and 75% with Playwright seeking, and check each story beat reads clearly.
- [ ] **Step 6: Commit** with the message "Studio: transcription and Meiosis I nondisjunction examples".

---

### Task 12: Documentation

**Files:**
- Modify:
  - `web/help.html`: add a "Studio" section covering the script format, actor and verb tables (generated by hand from `schema.js`), examples, export and the "not to scale" guidance.
  - `README.md`: a short "Animation Studio" section.
  - `CHANGELOG.md`: an "Unreleased" entry.
- Test: `tests/test_studio_static.py`

- [ ] **Step 1: Write the failing test.** `test_help_documents_every_verb`: parse the `VERBS` keys out of `web/studio/schema.js` with a regex over the object literal, and assert each appears in `help.html` inside a `<code>` element.
- [ ] **Step 2: Run it.** Expected: FAIL.
- [ ] **Step 3: Write the docs.**
- [ ] **Step 4: Run the full suite.** Run `python -m unittest discover -s tests`, `node --test web/studio/test/` and `python -m unittest e2e.test_smoke e2e.test_studio`. Expected: all PASS.
- [ ] **Step 5: Commit** with the message "Studio: help page, README and changelog", then push the branch.
