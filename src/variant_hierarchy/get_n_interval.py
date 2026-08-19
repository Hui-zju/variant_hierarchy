from __future__ import annotations

import sys
from functools import lru_cache

from .get_accession import DEFAULT_UTA_DSN, GenomicBackendUnavailable, get_accession_info
from .normalize.aa_position import parse_position_range


def _load_hgvs_modules():
    try:
        import hgvs.dataproviders.uta
        import hgvs.parser
        import hgvs.variantmapper
    except ModuleNotFoundError as exc:
        raise GenomicBackendUnavailable(
            "contains_variant requires hgvs and a local UTA service. Install the genomic/research extra."
        ) from exc
    return hgvs.dataproviders.uta, hgvs.parser, hgvs.variantmapper


@lru_cache(maxsize=8)
def _hgvs_context(uta_dsn: str):
    uta, parser_mod, variantmapper_mod = _load_hgvs_modules()
    try:
        hdp = uta.connect(uta_dsn)
    except Exception as exc:
        raise GenomicBackendUnavailable(
            "Unable to connect to UTA. contains_variant requires a local UTA database. "
            f"UTA DSN: {uta_dsn!r}."
        ) from exc
    hp = parser_mod.Parser()
    vm = variantmapper_mod.VariantMapper(hdp)
    return hdp, hp, vm


@lru_cache(maxsize=None)
def get_exons(gene_symbol: str, uta_dsn: str = DEFAULT_UTA_DSN):
    ac_info = get_accession_info(gene_symbol, uta_dsn=uta_dsn)
    if ac_info is None:
        return None
    tx_ac, alt_ac, aln = ac_info
    hdp, _, _ = _hgvs_context(uta_dsn)
    exons = hdp.get_tx_exons(tx_ac, alt_ac, aln)
    if not exons:
        raise ValueError(f"No exons returned for {gene_symbol} ({tx_ac}, {alt_ac}, {aln})")
    return tx_ac, exons


@lru_cache(maxsize=None)
def gene_to_n_interval(gene_symbol: str, uta_dsn: str = DEFAULT_UTA_DSN):
    exon_info = get_exons(gene_symbol, uta_dsn=uta_dsn)
    if not exon_info:
        return gene_symbol, 0, sys.maxsize
    tx_ac, exons = exon_info
    start = min(int(e["tx_start_i"]) for e in exons)
    end = max(int(e["tx_end_i"]) for e in exons)
    return tx_ac, start, end


@lru_cache(maxsize=None)
def exon_to_n_interval(gene_symbol: str, exon_pos: str, uta_dsn: str = DEFAULT_UTA_DSN):
    exon_info = get_exons(gene_symbol, uta_dsn=uta_dsn)
    if not exon_info:
        return None
    tx_ac, exons = exon_info
    a, b = parse_position_range(exon_pos)
    if a is None or b is None:
        return None

    ord_min = min(a, b) - 1
    ord_max = max(a, b) - 1
    chosen = [e for e in exons if ord_min <= int(e["ord"]) <= ord_max]
    if not chosen:
        return None

    start = min(int(e["tx_start_i"]) for e in chosen)
    end = max(int(e["tx_end_i"]) for e in chosen)
    return tx_ac, start, end


@lru_cache(maxsize=None)
def p_pos_to_n_interval(gene_symbol: str, p_pos: str, uta_dsn: str = DEFAULT_UTA_DSN):
    ac_info = get_accession_info(gene_symbol, uta_dsn=uta_dsn)
    if ac_info is None:
        return None
    tx_ac, alt_ac, aln = ac_info
    hdp, _, _ = _hgvs_context(uta_dsn)
    tx_info = hdp.get_tx_info(tx_ac, alt_ac, aln)
    cds_start_1based = tx_info["cds_start_i"] + 1
    cds_end_1based = tx_info["cds_end_i"]

    aa1, aa2 = parse_position_range(p_pos)
    if aa1 is None or aa2 is None:
        return None
    aa1, aa2 = sorted([aa1, aa2])
    if aa1 < 1:
        raise ValueError(f"Protein positions must be >=1, got {p_pos}")

    c_start = 3 * (aa1 - 1) + 1
    c_end = 3 * aa2
    n_start = cds_start_1based + c_start - 1
    n_end = min(cds_start_1based + c_end - 1, cds_end_1based)
    return tx_ac, n_start, n_end


@lru_cache(maxsize=None)
def p_hgvs_to_n_interval(gene_symbol: str, p_hgvs: str, uta_dsn: str = DEFAULT_UTA_DSN):
    ac_info = get_accession_info(gene_symbol, uta_dsn=uta_dsn)
    if ac_info is None:
        return None
    tx_ac, alt_ac, aln = ac_info
    hdp, hp, _ = _hgvs_context(uta_dsn)
    tx_info = hdp.get_tx_info(tx_ac, alt_ac, aln)
    cds_start_1based = tx_info["cds_start_i"] + 1
    cds_end_1based = tx_info["cds_end_i"]

    var_p = hp.parse_hgvs_variant(f"{tx_ac}:p.{str(p_hgvs).strip()}")
    pos = var_p.posedit.pos
    aa1 = getattr(pos.start, "pos", None)
    aa2 = getattr(pos.end, "pos", None)
    if aa1 is None or aa2 is None:
        raise ValueError(f"Unable to extract amino acid range from {p_hgvs!r}")

    aa1, aa2 = sorted([int(aa1), int(aa2)])
    if aa1 < 1:
        raise ValueError(f"Protein positions must be >=1, got {p_hgvs}")

    c_start = 3 * (aa1 - 1) + 1
    c_end = 3 * aa2
    n_start = cds_start_1based + c_start - 1
    n_end = min(cds_start_1based + c_end - 1, cds_end_1based)
    return tx_ac, n_start, n_end


@lru_cache(maxsize=None)
def c_hgvs_to_n_interval(gene_symbol: str, c_hgvs: str, uta_dsn: str = DEFAULT_UTA_DSN):
    ac_info = get_accession_info(gene_symbol, uta_dsn=uta_dsn)
    if ac_info is None:
        return None
    tx_ac, _, _ = ac_info
    _, hp, vm = _hgvs_context(uta_dsn)

    var_c = hp.parse_hgvs_variant(f"{tx_ac}:{c_hgvs}")
    var_n = vm.c_to_n(var_c)
    pos = var_n.posedit.pos
    start = pos.start.base + getattr(pos.start, "offset", 0)
    end = pos.end.base + getattr(pos.end, "offset", 0)
    return tx_ac, start, end


def to_n_interval(pos, uta_dsn: str = DEFAULT_UTA_DSN):
    gene = pos.get("gene")
    gene2 = pos.get("gene2")
    exon = pos.get("exon")
    amino_acid = pos.get("amino_acid")
    p_hgvs = pos.get("p_hgvs_posedit")
    c_hgvs = pos.get("c_hgvs_desc")

    if gene and exon:
        return [exon_to_n_interval(gene, exon, uta_dsn=uta_dsn)]
    if gene and amino_acid:
        return [p_pos_to_n_interval(gene, amino_acid, uta_dsn=uta_dsn)]
    if gene and p_hgvs:
        return [p_hgvs_to_n_interval(gene, p_hgvs, uta_dsn=uta_dsn)]
    if gene and c_hgvs:
        return [c_hgvs_to_n_interval(gene, c_hgvs, uta_dsn=uta_dsn)]
    if gene and gene2:
        return [
            gene_to_n_interval(gene, uta_dsn=uta_dsn),
            gene_to_n_interval(gene2, uta_dsn=uta_dsn),
        ]
    if gene:
        return [gene_to_n_interval(gene, uta_dsn=uta_dsn)]
    return None
