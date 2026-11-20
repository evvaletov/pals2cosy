# Changelog

## [0.6.1] — 2026-11-20

### Fixed
- Package metadata: the author e-mail address was wrong; it is now `evv@msu.edu`.
  No code changes.

## [0.6.0] — 2026-11-19

Targets the PALS specification at commit `b9d5b355` (13 November 2026),
recorded as `pals2cosy.__pals_commit__` and printed by `--version`. 0.5.0 was
not published on PyPI, so PyPI users upgrade from 0.4.2.

### Fixed
- **Multipole field scaling.** PALS `BnN` are field derivatives in T/m^N
  (`By + i Bx = Σ (BnN + i BsN)(x + i y)^N / N!`), not pole-tip fields. COSY's
  `MQ`/`MH`/`MO`/`M5` now receive `BnN · r^N / N!` at the bore radius `r`.
  Earlier releases passed `BnN` through unchanged, so a quadrupole gradient
  of 1 T/m became a pole-tip field of 1 T (a 74× too strong quadrupole at the
  default aperture). Lattices written for 0.5 with pole-tip values in `Bn1`
  must be rewritten with gradients.
- RF cavities: COSY's `RF` needs its voltage as an array and is a thin kick;
  the cavity is now emitted as `DL L/2 ; RF … ; DL L/2`, preserving the
  element length (earlier output dropped the length and failed in COSY 10.2).
- Bends with a negative angle: face angles are emitted in the `CB`-mirrored
  frame (`sign(θ)·e`); a rectangular negative bend previously received the
  wrong edge focusing.
- `carbon12` used the neutral-atom mass for the bare nucleus.
- **Wigglers.** COSY's `WI` takes the period λ as its second argument; earlier
  releases passed the wavenumber 2π/λ, so COSY modelled a nearly uniform
  dipole (horizontal dipole focusing, no vertical wiggler focusing, a path
  length longer than the element).
- Bends without `ApertureP` gave COSY a half-gap of 0: `--fr 3` stopped with
  a division by zero, and `--fr 1`/`2` silently produced a hard edge. With
  fringe fields on, the half-gap now follows from `edge1_int`/`edge2_int`
  (see Added) or, failing those, is the bore radius (`--quad-aperture / 2`),
  with a warning.
- The reference particle is always set with `RP E M Z` and the PALS (CODATA
  2018) mass, followed by `RPS 0 G` where the anomalous moment `G` is
  tabulated. The masses of COSY's presets `RPP`, `RPE`, `RPMU` and `RPPI`
  differ from it by 1e-7 (proton, electron) to 1.5e-6 (pion), so COSY's
  reference rigidity disagreed with the one normalized strengths were
  converted with (by 2.7e-6 in the IOTA one-turn map).
- `--twiss` uses a 2D map unless `--dim 3` is given, which needs RF; COSY's
  `GT` cannot handle the longitudinal plane of a 3D map without RF.
- Warnings about standalone flat-format wedges (`DPW`), truncated Enge
  coefficients and an ignored `--beamline` are `PALSConversionWarning`s, so
  `--strict` catches them.

### Added
- `kind: Bend` with BendP two-of-three resolution ({`g_ref`, `Bn0_ref`,
  `radius_ref`} × {`length`, `L_chord`, `L_rectangle`} × {`angle_ref`}),
  `e1`/`e2` and `e1_rect`/`e2_rect` per `ref_geometry`, face curvatures
  `h1`/`h2`, and `tilt_ref` as `RA` rotations.
- Bend fringe-field integrals `edge1_int`/`edge2_int`: with `--fr` 1–3 they
  set the half-gap of `DIL` to `edge_int / 0.3695`, which gives COSY's
  default Enge field the PALS integral, in preference to `ApertureP`. In the
  IOTA example the bends then give the design vertical tune to 0.003. At
  `--fr 0`, where COSY's bends have hard edges, they warn.
- Normalized (`KnN`, `KsN`), integrated (`BnNL`, `KnNL`) and tapered
  multipole strengths, skew components and `tiltN` (via COSY `MMS`), orders up
  to 10, and zero-length (thin) multipoles.
- Reference particle from the `BeginningEle`'s `ReferenceP` (`species_ref`,
  `E_tot_ref` or `pc_ref`); RF `dE_ref` updates the reference used for
  normalized strengths downstream.
- PALS/OpenPMD species parser (`pals2cosy.parse_species`): leptons, hadrons,
  pions, light nuclei, and ions such as `#3He+2` or `Au+79`.
- `RFP` (voltage or gradient·`L_active`, frequency or `harmonic_number`,
  phase with the `zero_phase` conventions, `dE_ref`, also without a
  voltage) and `SolenoidP` (`Bsol`, `Ksol`). A `dE_ref` that differs from
  the gain COSY gives the reference particle warns, since COSY's reference
  follows its own gain.
- Safe evaluation of PALS expressions in numeric fields: constants and
  variables, arithmetic, PALS constants and functions, element parameter
  references. Expressions that cannot be evaluated (a negative base with a
  fractional exponent, overflow, `factorial` beyond 170, nesting too deep for
  the parser) are one-line conversion errors.
- Negative `repeat` (reversed element order, per the specification without
  direction reversal); in-place Lattice branches; facility-level and
  parameter-group `inherit`, where a parameter may be modified in another
  form (`Bn1` for an inherited `Kn1`, `Bsol` for `Ksol`); `use` of the last
  Lattice by default.
- BeamLine `length`, also as a line-item override: a field-free drift fills
  the BeamLine up to its declared length (at the front when a negative repeat
  reverses it); a length shorter than the items raises
  `UnsupportedFeatureError`.
- `ReferenceChange` `dtime_ref` and `time_ref` enter the branch reference time
  that `harmonic_number` refers to; a `ReferenceChange` keeps its length.
- `Inf` and `-Inf` in expressions; `random()` and `random_gauss()` are drawn
  once per facility definition, so all instances of an element share a value;
  `nint` rounds halves away from zero; a constant or variable may be defined
  again with an equal value written differently (`2*pi`, `twopi`).
- `--solenoid-fringe`: the fringe-field length of COSY's soft-edge solenoid
  `CMS`, which defaults to the bore radius. At `--fr 0` a warning gives the
  resulting focusing deficit relative to the hard-edge solenoid (about D/L).
- `PALS.version` check (unknown major versions warn), `extension_labels`,
  `reminders`.
- Explicit handling of everything not converted: `UnsupportedFeatureError`
  for `set`, `post_expansion`, ABSOLUTE controllers, `include`, `load`,
  superposition, `placement`, direction reversal, offset `Patch`es and
  reference-energy changes; `PALSConversionWarning` for forks, extra
  branches, multipass lines, adjustments, RELATIVE controllers, girders,
  misalignments, electric multipoles, kicks, multipole fields on elements
  that do not convert them, unknown keys and kinds (including unknown
  components of the parameter groups pals2cosy reads), and elements without a
  readable strength. `include` is rejected at any level of the document, and
  extension blocks (`name: {extension: type, ...}`) are passed over at any
  level. A branch whose first element is not its `BeginningEle` warns.
- `--strict` (conversion warnings become errors) and `--version`; the CLI
  reports conversion errors as one line with exit status 1.
- Conformance suite: the PALS `examples/` tree at the pinned commit in
  `tests/pals_examples/` with a manifest of expected outcomes, run in COSY
  INFINITY at `FR 0` to `3`, plus tests that compare COSY transfer maps with
  the PALS element definitions. The COSY harness checks COSY's exit status,
  compilation failures and Fortran runtime errors and supplies `SYSCA.DAT`
  (`PALS2COSY_COSY_SYSCA`) for `FR 1` and `2`.

### Changed (breaking)
- Internal element dicts carry `bn`/`bs` (PALS coefficients in T/m^N) and
  `b_pole` (flat-format pole-tip fields) instead of `bn1`–`bn5`.
- An unknown reference species is an error instead of a proton-mass fallback.
- A `FloorShift` must have zero length, like a `Marker`.
- `convert()` takes `dimensions=None` (3, or 2 with `twiss`) and
  `solenoid_fringe`.
- Generated programs set the reference particle with `RP E M Z` and
  `RPS 0 G` instead of COSY's presets `RPP`, `RPE`, `RPMU` and `RPPI` (see
  Fixed).
- Every PALS solenoid warns at `--fr 0` unless `--solenoid-fringe` is at most
  1/1000 of its length, because COSY's `CMS` is soft-edged and focuses less
  than the hard-edge PALS solenoid; `--strict` therefore rejects solenoids
  by default.
- Without a `use` statement the last Lattice is converted (previously the
  last BeamLine), as the specification prescribes.
- `BeginningEle` and `Placeholder` no longer appear in the element list.

### Deprecated
- `kind: SBend`/`RBend` (read as `Bend`), `Lattice.particle`,
  `SolenoidP.Bz`, `RFCavityP`; each warns with `PALSDeprecationWarning`.

## [0.5.0] — 2026-05-02

### Added
- Numeric validation helpers (`_to_float`, `_to_nonneg_float`,
  `_to_positive_float`) applied to all PALS element parameters: lengths,
  `MagneticMultipoleP.Bn*`, `BendP.g_ref`/`e1`/`e2`, `WigglerP.peak_field`/
  `period`, `SolenoidP.Bz`, `RFCavityP.*`, `ApertureP.y_width`. Malformed
  numeric input now raises `ValueError` with the offending element/field name.
  NaN and infinity are rejected.
- `WigglerP.period > 0` guard in the converter (was `ZeroDivisionError`).
- `length_m` vs `s_end_m - s_start_m` consistency check in the flat parser.
- Particle mass database (`PARTICLE_MASSES_MEV` in `constants.py`) covering
  the same set of species the converter emits reference-particle commands for
  (electron, positron, proton, antiproton, muon ±, pion ±, deuteron, alpha,
  carbon12). `--particle` override is now case-insensitive for mass lookup.

### Changed (breaking)
- Negative `repeat` and the `direction` modifier in official PALS lattices
  now raise `ValueError` instead of warning. Reversed-direction expansion is
  not implemented; silent forward expansion was misleading.
- Overlapping elements in the flat parser now raise `ValueError` instead of
  warning. Producing physically impossible COSYScript was masking real input
  errors.
- Top-level access in both parsers now raises `ValueError` with a clear
  message instead of raw `KeyError`/`TypeError`.

### Documentation
- User-visible text now calls the language COSYScript or COSY INFINITY input
  throughout, replacing the informal name derived from the `.fox` extension.

## [0.4.0] — 2026-03-21

### Added
- **Sextupole support:** `MH L Bn2 R ;` command from `MagneticMultipoleP.Bn2`.
- **Solenoid support:** `CMS Bz D L ;` command from `SolenoidP.Bz`.
- **Expanded particle support:** muon, pion±, positron, antiproton, deuteron,
  alpha, carbon12 — auto-selects RPE/RPP/RPMU/RPPI or generic RP.
- **`--twiss` flag:** Opt-in GT Twiss extraction (requires periodic lattice).
- **FR 0/1/2/3 support:** Correct FC emission for FR 1 and 3 (Enge ODE modes);
  omitted for FR 2 (symplectic scaling).
- Public API exports: `convert`, `parse_lattice`, `parse_pals_lattice`.
- Zero-length element guards for MQ, MH, DIL, CMS.

### Fixed
- Shallow-copy bugs in BeamLine repeat expansion.
- Standalone dipole Enge symmetry (same coefficients on both faces).
- μ⁺/antimuon now uses generic RP with charge +1, not RPMU (which hardcodes −1).

### Changed
- Twiss extraction is no longer emitted by default (use `--twiss` to enable).
- Shared file loading extracted to `_io.py`.

## [0.3.0] — 2026-03-03

### Added
- **Official PALS format support:** `pals_parser.py` parses `PALS:` root with
  `facility:` definitions, `BeamLine` composition, `inherit:`, and `repeat:`.
- **Quadrupole Bn1 support:** `MagneticMultipoleP.Bn1` (pole-tip field in Tesla).
- **SBend with BendP:** `BendP.g_ref` (1/m), `e1`/`e2` (radians).
- **Auto-detection** from root key (`PALS:` vs `beamline:`).
- **`--beamline NAME`** to select a specific BeamLine from the facility.

## [0.2.0] — 2026-03-02

### Added
- Particle type support (`--particle electron|proton`).
- Element name comments in COSYScript output (`--no-comments` to disable).

## [0.1.0] — 2026-02-22

Initial release. PALS lattice parsing, dipole consolidation, COSYScript generation
with Enge fringe fields (FC/FD) and negative-angle CB wrapping.
