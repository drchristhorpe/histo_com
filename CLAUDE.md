# CLAUDE.md

Guidance for Claude Code when working in this repository.

## What this is

`histo_com` computes the mass-weighted centre of mass of 3D biological
structures (PDB/mmCIF). It's a Python library, a Click-based CLI
(`histo-com`), and a Claude skill wrapping that CLI (`skills/histo-com/`).
See [README.md](README.md) for user-facing usage and
[docs/PLAN.md](docs/PLAN.md) for the design rationale.

## Environment

- Python 3.14, managed with `uv`. Use `uv sync`, `uv run <cmd>`, `uv run pytest`.
- Don't invoke a bare `python`/`pip` — always go through `uv run` /
  `uv add` so the lockfile stays authoritative.

## Layout

```
src/histo_com/
  core.py        # load_structure(), centre_of_mass(), HistoCom, convenience wrappers
  selectors.py   # domain/residue selector string+iterable parsing (DomainRef, ResidueRef)
  cli.py         # Click CLI (entry point: histo-com)
tests/
  fixtures/      # real downloaded .cif/.pdb files used by the tests, keep committed
  test_core.py
  test_selectors.py
  test_cli.py
skills/histo-com/SKILL.md
```

## Key invariants — don't break these

- `HistoCom(path)` parses the structure **once** in `__init__`; every
  `com*` method must reuse `self.structure`/`self.model`. Never re-parse
  per domain/residue in a loop.
- `centre_of_mass()` requires its input iterable to be **homogeneous**
  (all Residues, all Chains, or all Atoms) — it walks Bio.PDB's
  disorder-unpacking (`get_unpacked_list`) in level order and mixing
  levels at the start breaks that traversal.
- Disordered atoms/residues (altlocs) are resolved via
  `get_unpacked_list()`, matching Biopython's own
  `Entity.center_of_mass()` exactly — including **both** altloc variants,
  mass-weighted but *not* occupancy-weighted. This is Biopython's own
  convention, not a bug; don't "fix" it to pick a single altloc without
  checking with the user, and don't write test fixtures/comparisons that
  naively iterate `for atom in residue` on structures with altlocs (see
  `7r7y_1_aligned.pdb`, which has real altlocs, e.g. chain A residue 2) —
  that skips disorder resolution and won't match.
- `--mode domains` selectors collapse each token to **one** combined COM;
  `--mode residues` selectors **expand** ranges to one COM per residue.
  This asymmetry is intentional (see README's selector grammar table) —
  don't "fix" one to match the other.
- Only the first model (`structure[0]`) is used everywhere.

## Testing

- `uv run pytest` — fixtures are already downloaded into
  `tests/fixtures/`; tests don't hit the network.
- Fixtures map directly to the use cases the tool was built against:
  `1hhk_1_abd.cif` (single chain, whole-structure), `1hhk_1_peptide.cif`
  (single chain, per-residue), `8gvi_1_aligned.cif` (multi-chain, domain
  selectors: `P`, `L:1-180`, `A,B`), `7r7y_1_aligned.pdb` (native PDB
  format, a genuine partial-range chain domain, and real altlocs). If you
  add behaviour, prefer extending these existing fixtures over adding new
  downloads.
- When changing selector parsing or COM math, add a case to
  `test_selectors.py`/`test_core.py` rather than only eyeballing CLI
  output.

## Scope

The CLI intentionally exposes exactly four options: `filename`, `--mode`,
`--domains`, `--residues`. Don't add output-format flags or other options
without checking with the user first — this was a deliberate constraint,
not an oversight.
