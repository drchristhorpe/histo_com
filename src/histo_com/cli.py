"""Command line interface for histo_com."""

from __future__ import annotations

import click

from histo_com.core import HistoCom, StructureError
from histo_com.selectors import (
    SelectorError,
    format_domain,
    format_residue,
    parse_domains,
    parse_residues,
)


def _format_vector(v) -> str:
    return ", ".join(f"{c:.4f}" for c in v)


@click.command()
@click.argument("filename", type=click.Path(exists=True, dir_okay=False))
@click.option(
    "--mode",
    "-m",
    type=click.Choice(["all", "domains", "residues"]),
    default="all",
    show_default=True,
    help="Compute one centre of mass for the whole structure, one per "
    "domain (chain or chain+range), or one per residue.",
)
@click.option(
    "--domains",
    "-d",
    default=None,
    help="Domain selector, required for --mode domains. "
    "Comma-separated chains and/or chain:range tokens, "
    "e.g. 'A', 'A,B', or 'L:1-180'.",
)
@click.option(
    "--residues",
    "-r",
    default=None,
    help="Residue selector, required for --mode residues. "
    "Comma-separated residue numbers and/or ranges, "
    "e.g. '5', '1-9', or 'A:1-9,B:12'.",
)
@click.option(
    "--output",
    "-o",
    type=click.Path(dir_okay=False, writable=True),
    default=None,
    help="Write a PDB file with a pseudo-atom (HETATM, resName COM) marking "
    "each computed centre of mass, for viewing alongside the structure.",
)
def main(
    filename: str,
    mode: str,
    domains: str | None,
    residues: str | None,
    output: str | None,
) -> None:
    """Compute the centre of mass of a 3D biological structure (PDB/mmCIF).

    FILENAME is the path to a .cif/.mmcif or .pdb/.ent structure file.
    """
    if mode == "domains" and not domains:
        raise click.UsageError("--domains is required when --mode domains")
    if mode == "residues" and not residues:
        raise click.UsageError("--residues is required when --mode residues")

    try:
        histo_com = HistoCom(filename)

        if mode == "all":
            com = histo_com.com()
            click.echo(f"all: {_format_vector(com)}")

        elif mode == "domains":
            refs = parse_domains(domains)
            coms = histo_com.com_by_domains(refs)
            for ref, com in zip(refs, coms):
                click.echo(f"{format_domain(ref)}: {_format_vector(com)}")

        else:
            refs = parse_residues(residues)
            coms = histo_com.com_by_residues(refs)
            for ref, com in zip(refs, coms):
                click.echo(f"residue {format_residue(ref)}: {_format_vector(com)}")

        if output:
            histo_com.write_com_pdb(output, mode=mode, domains=domains, residues=residues)
            click.echo(f"Wrote centre-of-mass marker(s) to {output}")

    except (StructureError, SelectorError) as exc:
        raise click.ClickException(str(exc)) from exc


if __name__ == "__main__":
    main()
