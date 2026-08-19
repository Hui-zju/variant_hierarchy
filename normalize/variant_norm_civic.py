import re
import json
import sys
import pandas as pd
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from variant_hierarchy.normalize.variant_parser import parse_variant
from source_mappings.loader import load_yaml_mapping
CIVIC_MAPPING = load_yaml_mapping("civic_variant_aliases.yaml")
CIVIC_VARIANT_TYPE_TO_CLASS = load_yaml_mapping("civic_variant_type_classes.yaml")


def split_civic_variant_types(variant_type):
    value = str(variant_type or "").strip()
    if not value or value == "N/A":
        return set()
    return {item.strip() for item in value.split(",") if item.strip()}


def infer_variant_class_from_civic(variant_type):
    classes = {
        CIVIC_VARIANT_TYPE_TO_CLASS[item]
        for item in split_civic_variant_types(variant_type)
        if item in CIVIC_VARIANT_TYPE_TO_CLASS
    }
    return classes.pop() if len(classes) == 1 else None


def norm_variant(variant, gene, variant_type):
    variant = CIVIC_MAPPING.get(variant, variant)
    variant = variant.replace("FS", 'fs').replace("DEL", 'del').replace("INS", 'ins').replace("DUP", 'dup')
    normalized_variant_class = infer_variant_class_from_civic(variant_type)

    # ----------------------------------------split variant------------------------------------------
    joiner = ' and '
    variant = variant.replace(' + ', joiner).replace('; ', joiner).strip()
    if joiner in variant:
        result = []
        for n in variant.split(joiner):
            result.extend(norm_variant(n.strip(), gene, variant_type))
        return result


    # V600E/K
    match = re.match(r"^([A-Za-z])(\d+)([A-Za-z/]+)$", variant)
    if match and "/" in variant:
        ref_aa = match.group(1)
        position = match.group(2)
        alt_aa_set = match.group(3)
        alt_aa_list = alt_aa_set.split("/")  # ? "/" ??
        elements = [f"{ref_aa}{position}{alt_aa}" for alt_aa in alt_aa_list]
        results = []
        for element in elements:
            results.extend(norm_variant(element.strip(), gene, variant_type))
        return results


    # ----------------------------------------norm variant------------------------------------------
    parsed_variant = parse_variant(variant, gene, variant_class=normalized_variant_class)
    if parsed_variant:
        return [parsed_variant]

    # raise Exception(f"Unable to parse variant from variantName (variantName={variant}, reference1={gene})")
    print(f"Unable to parse variant from variantName (variantName={variant}, reference1={gene})")
    return []




def get_variant():
    file_path = "../data/sources/CIVIC/01-Nov-2024-VariantSummaries.tsv"
    df = pd.read_csv(file_path, sep='\t')
    records = df.to_dict(orient='records')
    variants = []
    for record in records:
        items = norm_variant(str(record['variant']), str(record['feature_name']), str(record['variant_types']))
        for item in items:
            variants.append({**item, **{"source": ["CIVIC"]},
                             **{"raw": {"alteration": str(record['variant']), "gene": str(record['feature_name']), "type": str(record['variant_types'])}}})
    unique_variants = list({json.dumps(v, sort_keys=True): v for v in variants}.values())
    # gene_variants = [variant for variant in unique_variants if variant['variantClass'] == "GeneVariant"]
    # unresolved_variants = [variant for variant in unique_variants if variant['variantClass'] == "UnresolvedVariant"]
    with open("civic_variants.json", "w", encoding="utf-8") as f:  # ../data/variants/
        json.dump(unique_variants, f, ensure_ascii=False, indent=4)
    return unique_variants


if __name__ == '__main__':
    variants = get_variant()
    # variants = [item for variant in variants for item in variant]
    # dic = {'aa_change': 'R177ins', 'displayName': 'VHL R177ins', 'reference1': 'VHL'}
    # if dic in variants:
    #     print('error')

    # res = norm_variant('FGFR2-TACC3', 'FANCC') # Loss-of-function  FGFR2::AHCYL1
    # data = set([str(record['variant_types']) for record in records if str(record['last_review_date'])])
    # print('FGFR2-TACC3' in data)
    print('ok')
