"""Parsing of domain/residue selector strings and tokens.

See docs/PLAN.md (section 5) for the selector grammar.
"""

from __future__ import annotations

from dataclasses import dataclass


class SelectorError(ValueError):
    """Raised for malformed or ambiguous domain/residue selectors."""


@dataclass(frozen=True)
class ResidueRef:
    """A single residue, optionally scoped to a chain."""

    chain: str | None
    resseq: int


@dataclass(frozen=True)
class DomainPart:
    """A chain, or a residue range within a chain — one piece of a domain."""

    chain: str | None
    start: int | None = None
    end: int | None = None


@dataclass(frozen=True)
class DomainRef:
    """A domain: one or more parts (chains/ranges) whose combined mass
    defines a single centre of mass, e.g. a chain, a chain+range, or
    several chains treated as one unit (``"A,B"``)."""

    parts: tuple[DomainPart, ...]


def _split_chain(token: str) -> tuple[str | None, str]:
    """Split ``"A:1-180"`` into ``("A", "1-180")``; ``"1-180"`` into ``(None, "1-180")``."""
    if ":" in token:
        chain, _, rest = token.partition(":")
        chain = chain.strip()
        if not chain:
            raise SelectorError(f"Empty chain id in selector token {token!r}")
        return chain, rest.strip()
    return None, token.strip()


def _is_range(text: str) -> bool:
    if "-" not in text:
        return False
    # Distinguish "1-180" (a range) from a bare negative number, which
    # residue sequence numbers never legitimately are in these structures.
    left, _, right = text.partition("-")
    return left.strip().lstrip("+").isdigit() and right.strip().isdigit()


def _tokenize(spec: str) -> list[str]:
    tokens = [t.strip() for t in spec.split(",")]
    tokens = [t for t in tokens if t]
    if not tokens:
        raise SelectorError("Empty selector")
    return tokens


def _parse_domain_part(token) -> DomainPart:
    """Parse one domain *part*: a bare chain letter (``"A"``), a bare
    residue range (``"1-180"``), or a chain-scoped range (``"A:1-180"``) /
    chain-scoped single residue (``"A:12"``).
    """
    if isinstance(token, DomainPart):
        return token
    if isinstance(token, int):
        return DomainPart(chain=None, start=token, end=token)

    text = str(token).strip()
    if not text:
        raise SelectorError("Empty domain token")

    chain, rest = _split_chain(text)

    if not rest:
        if chain is None:
            raise SelectorError(f"Empty domain token {token!r}")
        return DomainPart(chain=chain)

    if _is_range(rest):
        start_s, _, end_s = rest.partition("-")
        start, end = int(start_s), int(end_s)
        if start > end:
            raise SelectorError(f"Invalid residue range {rest!r}: start > end")
        return DomainPart(chain=chain, start=start, end=end)

    if rest.lstrip("+").isdigit():
        resnum = int(rest)
        return DomainPart(chain=chain, start=resnum, end=resnum)

    if chain is None:
        # A bare alphabetic token with no ':' is a chain letter, e.g. "A".
        return DomainPart(chain=text)

    raise SelectorError(f"Could not parse domain token {token!r}")


def parse_domain(token) -> DomainRef:
    """Parse one domain into a :class:`DomainRef`.

    A domain is a single chain/range, or several comma-joined parts that
    are combined into one centre of mass, e.g. ``"A,B"`` (chains A and B
    together as one domain) or ``"A,B:1-50"``.
    """
    if isinstance(token, DomainRef):
        return token
    if isinstance(token, (int, DomainPart)):
        return DomainRef(parts=(_parse_domain_part(token),))

    text = str(token).strip()
    if not text:
        raise SelectorError("Empty domain")
    parts = tuple(_parse_domain_part(t) for t in _tokenize(text))
    return DomainRef(parts=parts)


def parse_domains(spec) -> list[DomainRef]:
    """Parse a domains selector.

    A **string** is a single domain: comma-separated parts are combined
    into one centre of mass (``"A,B"`` -> one domain spanning chains A and
    B). An **iterable** describes multiple, separate domains, one per
    element (``["A", "B"]`` -> two domains, each its own centre of mass);
    an element may itself be a compound string such as ``"A,B"``.
    """
    if isinstance(spec, str):
        return [parse_domain(spec)]
    tokens = list(spec)
    if not tokens:
        raise SelectorError("Empty selector")
    return [parse_domain(t) for t in tokens]


def parse_residue_token(token) -> list[ResidueRef]:
    """Parse one residue token into a list of :class:`ResidueRef` (ranges expand)."""
    if isinstance(token, ResidueRef):
        return [token]
    if isinstance(token, int):
        return [ResidueRef(chain=None, resseq=token)]
    if isinstance(token, tuple):
        if len(token) == 2:
            chain, resseq = token
            return [ResidueRef(chain=chain, resseq=int(resseq))]
        raise SelectorError(f"Invalid residue tuple {token!r}")

    text = str(token).strip()
    if not text:
        raise SelectorError("Empty residue token")

    chain, rest = _split_chain(text)

    if _is_range(rest):
        start_s, _, end_s = rest.partition("-")
        start, end = int(start_s), int(end_s)
        if start > end:
            raise SelectorError(f"Invalid residue range {rest!r}: start > end")
        return [ResidueRef(chain=chain, resseq=r) for r in range(start, end + 1)]

    if rest.lstrip("+").isdigit():
        return [ResidueRef(chain=chain, resseq=int(rest))]

    raise SelectorError(f"Could not parse residue token {token!r}")


def _format_domain_part(part: DomainPart) -> str:
    if part.start is None:
        return part.chain
    range_text = str(part.start) if part.start == part.end else f"{part.start}-{part.end}"
    return f"{part.chain}:{range_text}" if part.chain else range_text


def format_domain(ref: DomainRef) -> str:
    """Render a :class:`DomainRef` back to its canonical selector text."""
    return ",".join(_format_domain_part(p) for p in ref.parts)


def format_residue(ref: ResidueRef) -> str:
    """Render a :class:`ResidueRef` back to its canonical selector text."""
    return f"{ref.chain}:{ref.resseq}" if ref.chain else str(ref.resseq)


def parse_residues(spec) -> list[ResidueRef]:
    """Parse a residues selector: a comma-separated string or an iterable of tokens."""
    if isinstance(spec, str):
        tokens = _tokenize(spec)
    else:
        tokens = list(spec)
        if not tokens:
            raise SelectorError("Empty selector")
    refs: list[ResidueRef] = []
    for t in tokens:
        refs.extend(parse_residue_token(t))
    return refs
