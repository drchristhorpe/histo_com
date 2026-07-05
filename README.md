# histo_com

Compute the **centre of mass** of 3D biological structures (PDB/mmCIF files)
— for a whole structure, for one or more chains/domains, or per residue.

Built on [Biopython](https://biopython.org/), it ships as:

- a Python library — `import histo_com`
- a CLI tool — `histo-com`
- a [Claude Code / Claude Desktop skill](skills/histo-com/SKILL.md)

Requires Python 3.14+.

## Install

```bash
uv sync                 # dev environment, from a checkout
uv tool install .       # install the `histo-com` CLI globally
# or
pip install .
```

## CLI usage

```
histo-com FILENAME [--mode {all,domains,residues}] [--domains SPEC] [--residues SPEC] [--output PATH]
```

- `FILENAME` — a `.cif`/`.mmcif` or `.pdb`/`.ent` structure file.
- `--mode`, `-m` — `all` (default), `domains`, or `residues`.
- `--domains`, `-d` — selector string, required for `--mode domains`.
- `--residues`, `-r` — selector string, required for `--mode residues`.
- `--output`, `-o` — write a PDB file with a pseudo-atom marking each
  computed centre of mass (see [Marker PDB output](#marker-pdb-output)).

### Whole structure

```bash
$ histo-com 1hhk_1_abd.cif --mode all
all: -42.3648, 56.0308, 63.6705
```

### Per residue

`--residues` ranges expand to one centre of mass per residue:

```bash
$ histo-com 1hhk_1_peptide.cif --mode residues --residues 1-9
residue 1: -51.1112, 61.8689, 64.0741
residue 2: -47.0694, 61.3467, 62.3303
...
residue 9: -28.4389, 60.3628, 64.7192
```

### Per domain (chain, or chain + residue range)

`--domains` always computes **one** combined centre of mass per
invocation. A comma joins chains/ranges into that one domain — e.g. `A,B`
means "chains A and B together", not "chain A, then chain B":

```bash
$ histo-com 8gvi_1_aligned.cif --mode domains --domains P
P: -39.6074, 62.0183, 64.0741

$ histo-com 8gvi_1_aligned.cif --mode domains --domains L:1-180
L:1-180: -45.8994, 33.1512, 49.5114

$ histo-com 8gvi_1_aligned.cif --mode domains --domains A,B
A,B: -43.2339, 100.2589, 60.7451
```

To get results for several *separate* domains, run the command once per
domain (the Python library's `com_by_domains()` can also take an
iterable of separate domains in a single call — see
[Library usage](#library-usage)).

### Marker PDB output

Pass `--output`/`-o` to also write a small PDB file containing one
pseudo-atom per computed centre of mass — useful for viewing the result
alongside the structure in PyMOL, ChimeraX, etc. Each marker is a
`HETATM` residue named `COM`, on its own chain when the domain maps to
exactly one chain (chain `Z` for a multi-chain combined domain or the
whole-structure marker), or keeping the original residue numbering
(residues):

```bash
$ histo-com 8gvi_1_aligned.cif --mode domains --domains A,B --output com.pdb
A,B: -43.2339, 100.2589, 60.7451
Wrote centre-of-mass marker(s) to com.pdb

$ cat com.pdb
HETATM    1  COM COM Z   1     -43.234 100.259  60.745  1.00  0.00           C
```

## Selector grammar

Used by both `--domains`/`--residues` and their library equivalents. Each
selector is built from tokens of the form:

| Token | Meaning |
|---|---|
| `A` | chain `A` |
| `12` | residue 12 on the *default chain* (only valid if the structure has exactly one chain) |
| `1-180` | residues 1–180 on the default chain |
| `A:12` | residue 12 on chain `A` |
| `A:1-180` | residues 1–180 on chain `A` |

The two modes combine tokens differently:

- **`--mode domains`**: a comma-separated **string** is **one domain** —
  its tokens are combined into a single centre of mass (`A,B` = chains A
  and B together). To compute several *separate* domains, either run the
  command once per domain, or (library only) pass an iterable where each
  element is its own domain, e.g. `["A", "B"]` -> two separate results.
- **`--mode residues`**: tokens/ranges always **expand** — one centre of
  mass per residue, regardless of whether the input is a string or an
  iterable.

## Library usage

`HistoCom` parses a structure file once; every method below reuses that
same in-memory structure, so passing an iterator of chains or residues
only costs one parse.

```python
from histo_com import HistoCom

h = HistoCom("8gvi_1_aligned.cif")

h.com()                          # -> np.ndarray, whole-structure centre of mass
h.com_by_domains("P")            # -> [np.ndarray], one domain (chain P)
h.com_by_domains("A,B")          # -> [np.ndarray], ONE combined domain (chains A+B together)
h.com_by_domains(["A", "B"])     # -> [np.ndarray, np.ndarray], two SEPARATE domains

h.com_by_residues(range(1, 10))  # -> 9 np.ndarrays, one per residue
h.com_by_residues("A:1-9,B:12")  # selector strings also work

h.write_com_pdb("com.pdb", mode="domains", domains=["P", "L:1-180", "A", "B"])
```

Domain/residue arguments accept a selector string, an iterable of chain
letters, an iterable of residue numbers, or `DomainRef`/`ResidueRef`
objects from `histo_com.selectors`. For domains, remember: a **string**
is one (possibly multi-chain) domain; an **iterable** is several separate
domains, one result per element.

One-shot convenience functions are also available, each still parsing the
file only once internally:

```python
from histo_com import com_all, com_domains, com_residues

com_all("structure.cif")
com_domains("structure.cif", ["A", "B"])
com_residues("structure.cif", range(1, 10))
```

## Notes and limitations

- Centre of mass is **mass-weighted** (by element), not a geometric
  centroid, using Biopython's per-atom `mass` and its disorder-aware
  (altloc) unpacking logic.
- Only the **first model** in a file is used — sufficient for X-ray/cryo-EM
  structures; NMR ensembles are not averaged across models.
- A chain-less residue/domain token (e.g. `1-180` with no `A:` prefix) is
  only accepted for single-chain structures; otherwise an error asks you
  to specify the chain.

## Development

```bash
uv sync
uv run pytest
```

Test fixtures under `tests/fixtures/` are real structure files (both
mmCIF and PDB format) downloaded from
[coordinates.histo.fyi](https://coordinates.histo.fyi/).

See [docs/PLAN.md](docs/PLAN.md) for the design rationale and
[CHANGELOG.md](CHANGELOG.md) for release history.
