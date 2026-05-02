# Changelog

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
- Replaced "FOX" with "COSYScript" or "COSY INFINITY input" throughout
  user-visible text. COSYScript is the official COSY INFINITY domain-specific
  language; "FOX" is a `.fox`-extension shorthand, not an official name.

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
