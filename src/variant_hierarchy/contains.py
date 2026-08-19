from typing import Any, Optional

from .get_accession import DEFAULT_UTA_DSN, GenomicBackendUnavailable
from .match_alteration import alteration_match
from .match_function import function_match
from .match_location import location_match
from .normalize.variant_parser import parse_variant


def _as_variant_dict(
    variant: dict[str, Any] | str,
    reference: Optional[str] = None,
    parse_kwargs: Optional[dict[str, Any]] = None,
) -> Optional[dict[str, Any]]:
    if isinstance(variant, dict):
        return variant
    if isinstance(variant, str):
        return parse_variant(variant, reference=reference, **(parse_kwargs or {}))
    raise TypeError("variant must be a parsed variant dict or a variant string")


def contains_variant(
    parent: dict[str, Any] | str,
    child: dict[str, Any] | str,
    *,
    parent_reference: Optional[str] = None,
    child_reference: Optional[str] = None,
    parse_kwargs: Optional[dict[str, Any]] = None,
    uta_dsn: str = DEFAULT_UTA_DSN,
) -> bool:
    """
    Return whether ``parent`` contains ``child``.

    This is the full containment algorithm. It resolves gene/exon/amino-acid
    locations through UTA and therefore requires the genomic optional
    dependencies plus a reachable local UTA database.
    """
    parent_variant = _as_variant_dict(parent, parent_reference, parse_kwargs)
    child_variant = _as_variant_dict(child, child_reference, parse_kwargs)

    if not parent_variant or not child_variant:
        return False

    rel_alteration = alteration_match(parent_variant["alteration"], child_variant["alteration"])
    if rel_alteration == "unrelated":
        return False

    rel_function = function_match(parent_variant["function"], child_variant["function"])
    if rel_function == "unrelated":
        return False

    rel_location = location_match(parent_variant["location"], child_variant["location"], uta_dsn=uta_dsn)
    if rel_location == "unrelated":
        return False

    return rel_location == "contains" or rel_alteration == "contains" or rel_function == "contains"


__all__ = ["GenomicBackendUnavailable", "contains_variant"]
