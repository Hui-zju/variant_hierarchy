import re
from functools import lru_cache
from pathlib import Path
from typing import Optional, Tuple


DEFAULT_UTA_DSN = "postgresql://uta_admin:uta@localhost:5432/uta/uta_20241220"
DATA_DIR = Path(__file__).resolve().parent / "data"
MANE_DATA_PATH = DATA_DIR / "MANE.GRCh38.v1.5.summary.txt"
CHROM_ALIAS_PATH = DATA_DIR / "hg38.chromAlias.txt"


class GenomicBackendUnavailable(RuntimeError):
    pass


def _missing_dependency_error(package: str) -> GenomicBackendUnavailable:
    return GenomicBackendUnavailable(
        "contains_variant requires the genomic optional dependencies and a local UTA service. "
        f"Missing or unavailable dependency: {package!r}. Install with the genomic/research extra "
        "and pass uta_dsn if your UTA database is not at the default local URL."
    )


def _load_hgvs_uta():
    try:
        import hgvs.dataproviders.uta
    except ModuleNotFoundError as exc:
        raise _missing_dependency_error("hgvs") from exc
    return hgvs.dataproviders.uta


def _load_mane_transcript_mappings():
    try:
        from cool_seq_tool.sources import ManeTranscriptMappings
    except ModuleNotFoundError as exc:
        raise _missing_dependency_error("cool-seq-tool") from exc
    return ManeTranscriptMappings


def gene_to_mane_record(gene_symbol: str):
    ManeTranscriptMappings = _load_mane_transcript_mappings()
    mane_client = ManeTranscriptMappings(mane_data_path=MANE_DATA_PATH, from_local=True)
    mane_records = mane_client.get_gene_mane_data(gene_symbol)
    if not mane_records:
        return None

    return next(
        (r for r in mane_records if str(r.get("MANE_status", "")).strip() == "MANE Select"),
        mane_records[0],
    )


def _load_hg38_nc_set() -> set[str]:
    nc_re = re.compile(r"\bNC_\d+\.\d+\b")
    nc_set: set[str] = set()
    with CHROM_ALIAS_PATH.open(encoding="utf-8") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            nc_set.update(nc_re.findall(line))
    return nc_set


def _nc_version(nc: str) -> int:
    match = re.search(r"\.(\d+)$", nc)
    return int(match.group(1)) if match else 0


@lru_cache(maxsize=8)
def _uta_provider(uta_dsn: str):
    uta = _load_hgvs_uta()
    try:
        return uta.connect(uta_dsn)
    except Exception as exc:
        raise GenomicBackendUnavailable(
            "Unable to connect to UTA. contains_variant requires a local UTA database. "
            f"UTA DSN: {uta_dsn!r}."
        ) from exc


def transcript_to_genomic_accession(
    transcript_accession: str,
    uta_dsn: str = DEFAULT_UTA_DSN,
) -> Tuple[str, str]:
    uta_provider = _uta_provider(uta_dsn)
    mapping_options = uta_provider.get_tx_mapping_options(transcript_accession)
    if not mapping_options:
        raise ValueError(f"No UTA mapping options found for transcript_ac={transcript_accession!r}")

    hg38_nc = _load_hg38_nc_set()
    nc_opts = [o for o in mapping_options if str(o.get("alt_ac", "")).startswith("NC_")]
    best = None
    if nc_opts:
        in_hg38 = [o for o in nc_opts if str(o.get("alt_ac")) in hg38_nc]
        pool = in_hg38 or nc_opts
        splign_pool = [o for o in pool if o.get("alt_aln_method") == "splign"] or pool
        best = max(splign_pool, key=lambda o: _nc_version(str(o.get("alt_ac", ""))))

    mapping_option_record = best or mapping_options[0]
    genomic_ac = mapping_option_record.get("alt_ac")
    alignment_method = mapping_option_record.get("alt_aln_method")
    if not genomic_ac or not alignment_method:
        raise ValueError(f"Chosen mapping option missing alt_ac/alt_aln_method: {mapping_option_record}")
    return genomic_ac, alignment_method


@lru_cache(maxsize=None)
def get_accession_info(
    gene_symbol: str,
    uta_dsn: str = DEFAULT_UTA_DSN,
) -> Optional[Tuple[str, str, str]]:
    mane_record = gene_to_mane_record(gene_symbol)
    if not mane_record:
        return None

    tx_ac = mane_record["RefSeq_nuc"]
    alt_ac, aln = transcript_to_genomic_accession(tx_ac, uta_dsn=uta_dsn)
    return tx_ac, alt_ac, aln
