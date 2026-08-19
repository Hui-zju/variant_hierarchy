import csv
import re
from functools import lru_cache
from pathlib import Path
from typing import Optional


AA_POSITION_ALIAS_PATH = Path(__file__).resolve().parent.parent / "data" / "aa_position_aliases.tsv"
DEFAULT_SEQCAT_BASE_URL = "https://mtb.bioinf.med.uni-goettingen.de"


def parse_position_range(x):
    nums = [int(n) for n in re.findall(r"\d+", str(x))]
    if not nums:
        return None, None
    if len(nums) == 1:
        a = nums[0]
        return (a, a) if a >= 1 else (None, None)
    a, b = nums[0], nums[1]
    if a < 1 or b < 1:
        return None, None
    return (a, b) if a <= b else (a, a)


@lru_cache(maxsize=1)
def _load_position_aliases() -> dict[tuple[str, str], str]:
    aliases: dict[tuple[str, str], str] = {}
    with AA_POSITION_ALIAS_PATH.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            gene = (row.get("gene") or "").strip().upper()
            position = (row.get("position") or "").strip()
            reference_aa = (row.get("reference_aa") or "").strip().upper()
            if gene and position and reference_aa:
                aliases[(gene, position)] = f"{reference_aa}{position}"
    return aliases


def get_hgvs_p_aa_range(
    gene: str,
    pos: str,
    timeout: int = 60,
) -> Optional[str]:
    if not gene or not pos:
        return None

    pos_s = str(pos).strip()
    local_value = _load_position_aliases().get((str(gene).upper(), pos_s))
    if local_value:
        return local_value

    return _get_hgvs_p_aa_range_seqcat(gene, pos_s, timeout=timeout)


def _get_hgvs_p_aa_range_seqcat(
    gene: str,
    pos: str,
    timeout: int = 60,
) -> Optional[str]:
    a, b = parse_position_range(pos)
    if a is None or b is None:
        return None
    if a != b:
        return f"{a}_{b}"

    try:
        import requests
    except ModuleNotFoundError:
        return None

    url = f"{DEFAULT_SEQCAT_BASE_URL}/CCS/v1/get_protein_sequence/{gene}"
    try:
        response = requests.get(url, timeout=timeout)
        response.raise_for_status()
        payload = response.json()
        aa = payload[0][gene]["data"]["protein_sequence"][a - 1]
    except Exception:
        return None
    return f"{aa}{a}" if aa else None
