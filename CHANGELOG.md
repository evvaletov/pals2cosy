# Changelog

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
- Element name comments in FOX output (`--no-comments` to disable).

## [0.1.0] — 2026-02-22

Initial release. PALS lattice parsing, dipole consolidation, FOX code generation
with Enge fringe fields (FC/FD) and negative-angle CB wrapping.
