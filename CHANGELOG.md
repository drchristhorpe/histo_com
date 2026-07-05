# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

### Changed

- **Breaking**: `--mode domains` / `com_by_domains()` now treats a
  comma-separated **string** (e.g. `"A,B"`) as a single domain whose
  chains/ranges are combined into one centre of mass, rather than as
  multiple separate domains. To get several separate domains' worth of
  results in one call, pass an **iterable** instead (e.g. `["A", "B"]`),
  which still returns one COM per element; from the CLI, run the command
  once per domain. `DomainRef` is now a container of one or more
  `DomainPart`s rather than a single chain/range.

### Added

- `HistoCom.write_com_pdb()` / CLI `--output`, `-o`: writes a PDB file
  with one `HETATM` pseudo-atom (resName `COM`) per computed centre of
  mass, for viewing alongside the structure in a molecular viewer.
  Domain markers are placed on the domain's own chain; residue markers
  keep the original residue numbering; the whole-structure marker uses
  chain `Z`.
- `histo_com.selectors.format_domain()` / `format_residue()`: render a
  `DomainRef`/`ResidueRef` back to its canonical selector text, shared
  by the CLI's text output and the new marker labelling.

## [0.1.0] - 2026-07-05

### Added

- `histo_com` Python library: `HistoCom` class that parses a PDB/mmCIF
  structure once and computes centre of mass for the whole structure
  (`.com()`), per domain/chain (`.com_by_domains()`), or per residue
  (`.com_by_residues()`).
- `centre_of_mass()` core function: mass-weighted, disorder/altloc-aware
  (matches Biopython's own `Entity.center_of_mass()` exactly), built on
  Biopython's per-atom `mass` values.
- Selector grammar (`histo_com.selectors`) supporting chain letters,
  residue ranges, and chain-scoped ranges (e.g. `A`, `1-180`, `A:1-180`),
  accepted as either comma-separated strings or Python iterables.
- `histo-com` CLI (Click-based) with `--mode {all,domains,residues}`,
  `--domains`, and `--residues` options.
- Claude Code / Claude Desktop skill (`skills/histo-com/`) wrapping the CLI.
- Test suite (pytest) against real downloaded structures from
  coordinates.histo.fyi (`1hhk_1_abd.cif`, `1hhk_1_peptide.cif`,
  `8gvi_1_aligned.cif`, `7r7y_1_aligned.pdb` — a native PDB-format file
  used to verify a genuine partial-range chain domain).
- `README.md`, `CLAUDE.md`, and design plan (`docs/PLAN.md`).
