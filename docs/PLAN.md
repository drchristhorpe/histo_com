# histo_com — Design & Implementation Plan

## 1. Purpose

`histo_com` computes the (mass-weighted) **centre of mass** for 3D biological
structures held in PDB or mmCIF files. It ships as:

1. A Python library (`import histo_com`)
2. A CLI tool (`histo-com`), built with Click
3. A Claude Code / Claude Desktop skill that wraps the CLI

## 2. Tooling

- Python **3.14**, managed with **uv** (`uv venv`, `uv sync`, `uv run`, `uv build`)
- **Biopython** for structure parsing (`Bio.PDB.MMCIFParser` / `PDBParser`)
- **Click** for the CLI
- `pyproject.toml` with a `src/` layout, `uv_build` build backend, and a
  `[project.scripts]` entry point (`histo-com = "histo_com.cli:main"`)
- `pytest` for tests

## 3. Centre-of-mass algorithm

Biopython 1.87 assigns every `Atom` a `.mass` attribute at parse time
(derived from its element) and `Entity.center_of_mass()` already implements
mass-weighted averaging, correctly resolving disordered atoms/residues
(altlocs) via `get_unpacked_list()`.

We reuse that exact approach but generalise it to work over an arbitrary
*homogeneous* starting collection of entities (e.g. a filtered list of
`Residue` objects for a chain+range domain), not just a single entity's
children. This lets one routine serve whole-structure, whole-chain,
chain+range, and single-residue calculations uniformly:

```python
def centre_of_mass(entities: Iterable) -> np.ndarray:
    """Mass-weighted centre of mass over a homogeneous iterable of
    Bio.PDB entities (Residues, Chains, ...) or Atoms."""
```

The starting entities must all be at the same Bio.PDB "level" (all
Residues, all Chains, ...) so the disorder-unpacking traversal stays
correct (this always holds for our call sites).

## 4. Structure loading

`histo_com.core.load_structure(path)`:
- picks `MMCIFParser` or `PDBParser` based on file extension (`.cif`/`.mmcif`
  vs `.pdb`/`.ent`), case-insensitively
- `QUIET=True` (suppress Biopython's construction warnings by default)
- always operates on the **first model** (`structure[0]`) — sufficient for
  X-ray/cryo-EM structures, documented as a limitation for NMR ensembles

## 5. Selector grammar (domains / residues)

Both CLI options and the library accept either a comma-separated string
(CLI-friendly) or a Python iterable (library-friendly) of tokens.

A token is one of:

| Form | Meaning |
|---|---|
| `A` | chain letter `A` |
| `1-180` | residue range 1–180 on the **default chain** (only valid if the structure has exactly one chain) |
| `A:1-180` | residue range 1–180 on chain `A` |
| `12` | single residue 12 on the default chain |
| `A:12` | single residue 12 on chain `A` |

- **`--mode domains`**: each token yields exactly **one** combined COM
  (e.g. `A:1-180` → one COM over that whole range; `A,B` → two COMs, one
  per chain).
- **`--mode residues`**: ranges expand element-wise — `1-9` yields **nine**
  COMs, one per residue, in order.
- **`--mode all`**: no selector needed; one COM over every atom in the
  first model.

Errors are raised (not silently guessed) when: a bare, chain-less token is
used on a multi-chain structure; a referenced chain/residue doesn't exist;
or the required option for the chosen mode is missing.

## 6. Library API (`histo_com/core.py`)

```python
from histo_com import HistoCom

h = HistoCom("structure.cif")        # parses once

h.com()                               # -> np.ndarray, mode "all"
h.com_by_domains(["A", "B"])          # -> list[np.ndarray]
h.com_by_domains(["L:1-180"])         # -> list[np.ndarray] (len 1)
h.com_by_residues(range(1, 10))       # -> list[np.ndarray] (len 9)
```

`HistoCom` parses the file once in `__init__`; all methods reuse the same
in-memory structure, satisfying "only loading the structure once" when
callers pass an iterator of chains/residues.

Module-level convenience functions (`histo_com.com_all(path)`, etc.) are
provided as thin wrappers for one-shot scripting use, each still only
parsing the file once internally.

## 7. CLI (`histo_com/cli.py`)

```
histo-com FILENAME --mode {all,domains,residues} [--domains SPEC] [--residues SPEC]
```

- `FILENAME`: positional argument
- `--mode` / `-m`: `all` (default), `domains`, `residues`
- `--domains` / `-d`: selector string, required when `--mode domains`
- `--residues` / `-r`: selector string, required when `--mode residues`

Output: one line per result, `<label>: x, y, z` (Å, 4 dp). This keeps the
CLI to exactly the four requested options — no extra flags.

## 8. Claude skill

`skills/histo-com/SKILL.md` — a skill definition describing when/how to
invoke the `histo-com` CLI (computing centres of mass for chains, domains,
or residue ranges in a structure file), so Claude Code/Desktop can call it
via the Bash tool once the package is installed (`uv tool install .` or
`pip install histo_com`).

## 9. Package layout

```
histo_com/
  pyproject.toml
  README.md
  CLAUDE.md
  CHANGELOG.md
  docs/PLAN.md
  src/histo_com/
    __init__.py
    core.py        # centre_of_mass(), HistoCom, load_structure
    selectors.py    # token parsing for domains/residues specs
    cli.py          # Click CLI
  skills/histo-com/SKILL.md
  tests/
    fixtures/       # downloaded test .cif files
    test_core.py
    test_selectors.py
    test_cli.py
```

## 10. Testing plan

Unit tests with pytest against small, real downloaded fixtures:

- `1hhk_1_abd.cif` — `all` mode
- `1hhk_1_peptide.cif` — `all` mode; residues iterator `range(1, 10)`
- `8gvi_1_aligned.cif` — domains `"P"`, `"L:1-180"`, `"A,B"`

Both library-level and CLI-level (via `click.testing.CliRunner`) tests.
Verify: correct number of results returned per call, coordinates are
finite 3-vectors, single-parse guarantee (structure loaded once per
`HistoCom` instance, checked via parser call count).

## 11. Workflow

1. Write this plan, commit it.
2. Scaffold project with `uv init`/`uv add`.
3. Implement `core.py`, `selectors.py`, `cli.py`.
4. Download fixtures, write and run tests.
5. Write README.md, CLAUDE.md, CHANGELOG.md.
6. Manually exercise the three documented use cases end-to-end.
