---
name: histo-com
description: Compute the centre of mass of a 3D biological structure (PDB/mmCIF file) — for the whole structure, for one or more chains/domains, or per residue. Use when asked for a structure's centre of mass, geometric centre, or "COM", or to locate a chain/domain/residue in 3D space by its mass-weighted centre.
---

# histo-com

`histo-com` is a CLI tool (installed from the `histo_com` package) that
computes the mass-weighted centre of mass of a PDB or mmCIF structure
file using Biopython. Invoke it with the Bash tool.

## When to use this skill

The user provides (or references) a `.cif`/`.mmcif` or `.pdb`/`.ent`
structure file and asks for its centre of mass, or the centre of mass of
a specific chain, a chain plus a residue range, or one or more individual
residues.

## Checking availability

```bash
histo-com --help
```

If this fails with "command not found", install it first:

```bash
uv tool install histo_com   # or: pip install histo_com
```

(If working from a checkout of the `histo_com` source repo instead of an
installed package, use `uv run histo-com ...` there instead.)

## Usage

```
histo-com FILENAME [--mode {all,domains,residues}] [--domains SPEC] [--residues SPEC]
```

- **Whole structure** (default mode): `histo-com structure.cif`
- **One or more chains**: `histo-com structure.cif --mode domains --domains A,B`
- **A chain restricted to a residue range** (one combined COM):
  `histo-com structure.cif --mode domains --domains A:1-180`
- **Individual residues** (one COM per residue; ranges expand):
  `histo-com structure.cif --mode residues --residues 1-9`
  `histo-com structure.cif --mode residues --residues A:1-9,B:12`

Selector rule of thumb: `--domains` tokens each collapse to **one** COM
(good for "where is chain A" or "where is this domain"); `--residues`
ranges **expand** to one COM per residue (good for "COM at each position
along this range").

A bare residue number/range with no `chain:` prefix (e.g. `--residues
1-9`) only works if the structure has exactly one chain — otherwise
`histo-com` reports the available chains and asks for a `chain:` prefix.

Output is one line per result: `<label>: x, y, z` in Ångströms.

## Example

```bash
$ histo-com 8gvi_1_aligned.cif --mode domains --domains P
P: -39.6074, 62.0183, 64.0741

$ histo-com 8gvi_1_aligned.cif --mode residues --residues 1-3
residue 1: ...
residue 2: ...
residue 3: ...
```

Report the coordinates back to the user in whatever form they asked for
(raw numbers, a table, etc.) — this skill only tells you how to obtain
them.
