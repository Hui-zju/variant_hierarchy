import re
from typing import Optional
from dataclasses import asdict
from .alteration_norm import get_normalized_mechanism
from .alteration_norm import get_normalized_function
from .variant import build_p_hgvs_variant, build_amino_acid_variant, build_exon_variant, build_fusion_variant, build_gene_variant, build_unresolved_variant

from .gene_norm import gene_norm
from .hgvs_parse import is_aa_mut, get_aa_mut_type, is_cds_mut, get_cds_mut_type


def parse_amino_acid_hgvs_variant(
    variant: str,
    gene: str,
    variant_class: Optional[str] = None,
):
    variant = variant.strip()
    variant = re.sub(r"^[pP]\.", "", variant)
    variant = variant.replace("DELINS", "delins").replace("DEL", "del").replace("INS", "ins")
    variant = variant.replace("→", ">").replace("‐", "-").replace("–", "-")
    n_gene = gene_norm(gene)
    # F530delF -> F530del; L755_T759delLRENT -> L755_T759del; keep delins unchanged.
    modified_variant = re.sub(
        r'^([A-Za-z]\d+(?:_[A-Za-z]\d+)?)(del|dup)(?!ins)[A-Za-z*]+$',
        r'\1\2',
        variant,
        flags=re.IGNORECASE,
    )
    if is_aa_mut(modified_variant):
        variant = modified_variant

    # Protein + CDS notation: 'R200W (c.598C>T)' or 'R474C c.1420C>T'
    match = re.match(r"^(?:p\.)?([a-z*]\d+\S*)\s+(?:\(c\.[^)]+\)|c\.\S+)$", variant, re.IGNORECASE)
    if match:
        variant = match.group(1)

    if is_aa_mut(variant) and n_gene:
        var_type = get_normalized_mechanism(
            get_aa_mut_type(variant),
            location="Amino Acid",
            variant_class=variant_class,
        )
        if var_type:
            return build_p_hgvs_variant(n_gene, variant, var_type)
    return None


def parse_amino_acid_categorical_variant(
    variant: str,
    gene: str,
    variant_class: Optional[str] = None,
):
    variant = variant.strip().strip("()")

    n_gene = gene_norm(gene)
    # ESR1 538
    match = re.match(r"^\d+$", variant)
    if n_gene and match:
        aa_range = variant
        var_type = get_normalized_mechanism("mutation", location="Amino Acid", variant_class=variant_class)
        return build_amino_acid_variant(n_gene, aa_range, var_type)

    # KIT (550-592)
    match = re.match(r"^(\d+)\s*[-\u2013\u2014]\s*(\d+)$", variant)
    if n_gene and match:
        aa_range = f"{match.group(1)}-{match.group(2)}"
        var_type = get_normalized_mechanism("mutation", location="Amino Acid", variant_class=variant_class)
        return build_amino_acid_variant(n_gene, aa_range, var_type)

    # amino acid 178-200 mutation
    match = re.match(r"^amino acid (\d+)(?:[-\u2013\u2014](\d+))?(?:\s+(.+))?$", variant, re.IGNORECASE)
    if n_gene and match:
        pos1 = match.group(1)
        pos2 = match.group(2)
        pos = pos1 + "-" + pos2 if pos2 else pos1
        raw_var_type = match.group(3) or "mutation"
        var_type = get_normalized_mechanism(raw_var_type, location="Amino Acid", variant_class=variant_class)
        if var_type:
            return build_amino_acid_variant(n_gene, pos, var_type)

    # X1473_splice
    match = re.match(r"^(?:p\.)?[a-zA-Z*](\d+)_(.+)$", variant, re.IGNORECASE)
    if n_gene and match:
        pos = match.group(1)
        var_type = get_normalized_mechanism(match.group(2), location="Amino Acid", variant_class=variant_class)
        if var_type:
            return build_amino_acid_variant(n_gene, pos, var_type)

    # G12 mutation
    match = re.match(r"^[a-zA-Z](\d+)\s+(.+)$", variant, re.IGNORECASE)
    if n_gene and match:
        pos = match.group(1)
        var_type = get_normalized_mechanism(match.group(2), location="Amino Acid", variant_class=variant_class)
        if var_type:
            return build_amino_acid_variant(n_gene, pos, var_type)

    # T1151ins / P780ins: position-level shorthand, not complete p.HGVS.  # |del|dup
    match = re.match(r"^(?:p\.)?[a-zA-Z](\d+)(ins)$", variant, re.IGNORECASE)
    if n_gene and match:
        pos = match.group(1)
        var_type = get_normalized_mechanism(match.group(2), location="Amino Acid", variant_class=variant_class)
        if var_type:
            return build_amino_acid_variant(n_gene, pos, var_type)

    # V600   V600X   is_aa_mut('V600X') = True
    match = re.match(r"^(?:p\.)?[a-zA-Z](\d+)X?$", variant)
    if n_gene and match:
        aa_range = match.group(1)
        var_type = get_normalized_mechanism("missense mutation", location="Amino Acid", variant_class=variant_class)
        return build_amino_acid_variant(n_gene, aa_range, var_type)


    # 'PDGFRA 311_848del'
    match = re.match(r'(\d+)_(\d+)(del|ins)$', variant)
    if n_gene and match:
        pos1 = match.group(1)
        pos2 = match.group(2)
        aa_range = pos1 + "-" + pos2
        var_type = get_normalized_mechanism(match.group(3), location="Amino Acid", variant_class=variant_class)
        return build_amino_acid_variant(n_gene, aa_range, var_type)


    # ERBB2 inframe deletion (755-759)
    match = re.match(r'^(.+)?\s?\((\d+-\d+)\)$', variant)
    if n_gene and match:
        var_type = match.group(1)
        aa_range = match.group(2)
        var_type = var_type if var_type else "mutation"
        var_type = get_normalized_mechanism(var_type, location="Amino Acid", variant_class=variant_class)
        if var_type:
            return build_amino_acid_variant(n_gene, aa_range, var_type)

    # NOTCH1 2245_2536 splice acceptor variant
    match = re.match(r'^(\d+_\d+)\s(.+)$', variant)
    if n_gene and match:
        aa_range = match.group(1)
        var_type = match.group(2)
        var_type = var_type if var_type else "mutation"
        var_type = get_normalized_mechanism(var_type, location="Amino Acid", variant_class=variant_class)
        if var_type:
            return build_amino_acid_variant(n_gene, aa_range, var_type)

    return None



def parse_exon_variant(
    variant: str,
    gene: str,
    variant_class: Optional[str] = None,
):
    n_gene = gene_norm(gene)

    if not n_gene:
        return None

    variant = re.sub(r"\s*\((?:EGFR\s*)?v?III\)\s*$", "", variant, flags=re.IGNORECASE).strip()
    exon_pattern = re.compile(
        r"^exons?[\s_]*(\d+)"
        r"(?:[\s_]*[-\u2013\u2014][\s_]*(\d+))?"
        r"(?:[\s_]*([A-Za-z].*))?"
        r"\s*$",
        re.IGNORECASE,
    )

    match = exon_pattern.match(variant)
    if not match:
        return None

    exon_start = match.group(1)
    exon_end = match.group(2)
    exon = f"{exon_start}-{exon_end}" if exon_end else exon_start
    raw_var_type = (match.group(3) or "mutation").strip()

    var_type = get_normalized_mechanism(
        raw_var_type,
        location="Exon",
        variant_class=variant_class,
        gene=n_gene,
        exon=exon,
    )
    var_func = get_normalized_function(raw_var_type)

    if var_type or var_func:
        return build_exon_variant(n_gene, exon, var_type, var_func)
    return None


def parse_gene_variant(
    variant: str,
    gene: str,
    variant_class: Optional[str] = None,
):
    n_gene = gene_norm(gene)
    var_type = get_normalized_mechanism(variant, location="Gene", variant_class=variant_class)
    var_func = get_normalized_function(variant)
    if n_gene and (var_type or var_func):
        return build_gene_variant(n_gene, var_type, var_func)
    return None

def parse_fusion_variant(
    variant: str,
    gene: Optional[str] = None,
    variant_class: Optional[str] = None,
):
    hyphen_pattern = r"[\-–—_]"
    variant = variant.strip()
    n_gene = gene_norm(gene)
    var_type = get_normalized_mechanism(variant, location="Gene", variant_class=variant_class)

    # variant: fusion, gene: ALK
    if n_gene and var_type == "fusion":
        return build_fusion_variant(n_gene)


    # variant: CDKN2B-AS1-PTPRD fusion, gene: CDKN2B-AS1
    mutil_hyphen_match = re.match(r'^(?P<core>[^\s]+-[^\s]+)\s+fusions?$', variant, re.IGNORECASE)
    if n_gene and mutil_hyphen_match:
        core = mutil_hyphen_match.group("core").strip()
        if core.startswith(gene):
            g2 = core[len(gene):].strip("-")
            n_g2 = gene_norm(g2)
            if n_g2:
                return build_fusion_variant(n_gene, n_g2)
        if core.endswith(gene):
            g1 = core[:-len(gene)].strip("-")
            n_g1 = gene_norm(g1)
            if n_g1:
                return build_fusion_variant(n_g1, n_gene)


    # variant: ALK fusion ,  ALK-EGFR fusion
    hyphen_match = re.match(
        rf"^([A-Za-z0-9]+){hyphen_pattern}?([A-Za-z0-9]+)?{hyphen_pattern}?\s*fusions?$",
        variant,
        re.IGNORECASE
    )
    if hyphen_match:
        g1 = hyphen_match.group(1)
        g2 = hyphen_match.group(2)
        n_g1 = gene_norm(g1)
        n_g2 = gene_norm(g2)
        if n_g1 and n_g2:
            return build_fusion_variant(n_g1, n_g2)
        if n_gene and n_g1:
            if n_gene == n_g1 and g2 is None:
                return build_fusion_variant(n_g1)
            if n_gene != n_g1 and g2 is None:
                return build_fusion_variant(n_gene, n_g1)


    # variant: EGFR-ALK     (gene: CDKN2B-AS1)
    match = re.match(
        rf"^([a-z0-9]+){hyphen_pattern}([a-z0-9]+)$",
        variant,
        re.IGNORECASE
    )
    if match:
        g1, g2 = match.groups()
        n_g1 = gene_norm(g1)
        n_g2 = gene_norm(g2)
        if n_g1 and n_g2:
            return build_fusion_variant(n_g1, n_g2)

    # variant: IGH::BCL2  IGH::?   v::BCL2
    double_colon_match = re.match(
        r"^([a-z0-9\-?]+)::([a-z0-9\-?]+)(?:\s+fusions?)?$",
        variant,
        re.IGNORECASE
    )
    if double_colon_match:
        g1, g2 = double_colon_match.groups()
        n_g1 = gene_norm(g1)
        n_g2 = gene_norm(g2)
        if n_g1 and n_g2:
            return build_fusion_variant(n_g1, n_g2)
        if n_g1 and g2 in ["?", "v"]:
            return build_fusion_variant(n_g1, g2)
        if n_g2 and g1 in ["?", "v"]:
            return build_fusion_variant(g1, n_g2)

    # variant: fusion   gene:IGH::BCL2  EGFR-ALK
    if gene and var_type == "fusion":
        match = re.match(
            rf"^([a-z0-9]+){hyphen_pattern}([a-z0-9]+)$",
            gene,
            re.IGNORECASE
        )
        if match:
            g1, g2 = match.groups()
            n_g1 = gene_norm(g1)
            n_g2 = gene_norm(g2)
            if n_g1 and n_g2:
                return build_fusion_variant(n_g1, n_g2)

        double_colon_match = re.match(
            r"^([a-z0-9\-?]+)::([a-z0-9\-?]+)(?:\s+fusions?)?$",
            gene,
            re.IGNORECASE
        )
        if double_colon_match:
            g1, g2 = double_colon_match.groups()
            n_g1 = gene_norm(g1)
            n_g2 = gene_norm(g2)
            if n_g1 and n_g2:
                return build_fusion_variant(n_g1, n_g2)
            if n_g1 and g2 in ["?", "v"]:
                return build_fusion_variant(n_g1, g2)
            if n_g2 and g1 in ["?", "v"]:
                return build_fusion_variant(g1, n_g2)
    return None


def parse_unsolved_variant(variant: str, gene: Optional[str] = None):
    if gene:
        print(f"Unable to parse variant from variantName (variantName={variant}, reference1={gene})")
    else:
        print(f"Unable to parse variant from variantName (variantName={variant})")
    return build_unresolved_variant(variant, gene)


VARIANT_TEXT_ALIASES = {
    "EGFRVIII": "EGFR exon 2-7 deletion",
    "VIII": "exon 2-7 deletion",
}

def preprocess_variant_text(text: Optional[str], apply_alias: bool = False) -> Optional[str]:
    if text is None:
        return None

    text = text.strip()
    text = re.sub(r"\s+", " ", text)
    text = re.sub(r"\bmutations\b", "mutation", text, flags=re.IGNORECASE)
    text = re.sub(r"\bvariants\b", "variant", text, flags=re.IGNORECASE)
    text = re.sub(r"\bfusions\b", "fusion", text, flags=re.IGNORECASE)
    text = re.sub(r"\bdeletions\b", "deletion", text, flags=re.IGNORECASE)
    text = re.sub(r"\binsertions\b", "insertion", text, flags=re.IGNORECASE)
    text = re.sub(r"\bduplications\b", "duplication", text, flags=re.IGNORECASE)
    if apply_alias:
        alias_key = re.sub(r"[^A-Z0-9]", "", str(text or "").upper())
        text = VARIANT_TEXT_ALIASES.get(alias_key, text)
    return text


def parse_variant(
    variant: str,
    reference: Optional[str] = None,
    variant_class: Optional[str] = None,
):
    """
    General variant annotation entry point.
    - variant: Single variant description or gene+variant, e.g. ‘BRAF V600E’ / “V600E” / ‘BCR-ABL1’
    - gene: Passed if already extracted from variant; otherwise may be None
    """
    if not variant:
        return None

    variant = preprocess_variant_text(variant, apply_alias=True)
    reference = preprocess_variant_text(reference)

    if not reference:
        split = re.match(r'^([A-Za-z0-9.\-/]+)[_\s]+(.*)$', variant)
        if split:
            reference = split.group(1)
            variant = preprocess_variant_text(split.group(2), apply_alias=True)

    # Parsing order：fusion → AA change → AA region → position-only → exon → gene-level
    parsed = (
        parse_amino_acid_categorical_variant(
            variant,
            reference,
            variant_class,
        )
        or parse_amino_acid_hgvs_variant(
            variant,
            reference,
            variant_class,
        )
        or parse_exon_variant(
            variant,
            reference,
            variant_class,
        )
        or parse_fusion_variant(
            variant,
            reference,
            variant_class,
        )
        or parse_gene_variant(
            variant,
            reference,
            variant_class,
        )
        # or parse_unsolved_variant(variant, reference, variant_class)
    )
    if parsed:
        return asdict(parsed)
    else:
        return None

if __name__ == '__main__':
    # res= parse_fusion_variant("1p/19g)")
    # 1p/19q deletion   10p12.3 rearrangement  17p deletion
    # res= parse_chromosome_variant("deletion", "1p/19q")
    # RAD21 RAD21-C8orf37-AS1 fusion , PTPRD CDKN2B-AS1-PTPRD fusion

    # res= parse_variant("NF1 intragenic")
    # if res:
    #     res_dict = res.to_dict()
    # print('ok')
    # variant = 'Splice Site (c.3028G>A)'
    # match = re.match(r'\((.*?)\)', variant)
    # if match:
    #     variant = match.group(1)
    #     print('ok')
    # res = parse_nucleotide_hgvs_variant("Splice Site (c.3028G>A)", "EGFR")

    # res = parse_amino_acid_categorical_variant("K618del", "BRAF")
    res = parse_variant("EGFR exon 20 insertions")
    if res:
        print(res)
