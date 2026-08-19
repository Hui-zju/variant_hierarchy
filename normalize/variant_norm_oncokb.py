import re
import json
import sys
import pandas as pd
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from variant_hierarchy.normalize.variant_parser import parse_variant
from source_mappings.loader import load_yaml_mapping
ONCOKB_MAPPING = load_yaml_mapping("oncokb_alteration_aliases.yaml")

def norm_variant(variant, gene):
    variant = ONCOKB_MAPPING.get(variant.lower().strip(), variant).replace("p.", "")

    # ----------------------------------------split variant------------------------------------------
    # Oncogenic Mutations (excluding V600)    V600 (excluding V600E and V600K)
    pattern = r"^(.+)\s\(excluding\s([^\)]+)\)"
    matches = re.search(pattern, variant)
    if matches:
        element = matches.group(1)
        # Extract the list of mutations
        excluded_elements = matches.group(2).replace("and", ",").split(",")
        # Clean up extra spaces
        excluded_elements = [mutation.strip() for mutation in excluded_elements]
        all_elements = [element] + excluded_elements

        results = []
        for element in all_elements:
            results.extend(norm_variant(element.strip(), gene))
        return results

    # D538, E380, L469V, L536, S463P, Y537
    joiner = ','
    if joiner in variant:
        result = []
        for n in variant.split(joiner):
            result.extend(norm_variant(n.strip(), gene))
        return result

    # ----------------------------------------norm variant------------------------------------------
    parsed_variant = parse_variant(variant, gene, variant_class=None)
    if parsed_variant:
        return [parsed_variant]

    # raise Exception(f"Unable to parse variant from variantName (variantName={variant}, gene={gene})")
    print(f"Unable to parse variant from variantName (variantName={variant}, gene={gene})")
    return []

def get_variant():
    file_path = "../data/sources/Oncokb/oncokb_biomarker_drug_associations.tsv"
    df = pd.read_csv(file_path, sep='\t')
    records = df.to_dict(orient='records')
    variants = []
    for record in records:
        items = norm_variant(record['Alterations'], record['Gene'])
        for item in items:
            variants.append({**item, **{"source": ["OncoKB"]}, **{"raw": {"alteration": record['Alterations'], "gene": record['Gene']}}})
    unique_variants = list({json.dumps(v, sort_keys=True): v for v in variants}.values())
    # gene_variants = [variant for variant in unique_variants if variant['variantClass'] == "GeneVariant"]
    # unresolved_variants = [variant for variant in unique_variants if variant['variantClass'] == "UnresolvedVariant"]
    with open("oncokb_variants.json", "w", encoding="utf-8") as f:  # ../data/variants/
        json.dump(unique_variants, f, ensure_ascii=False, indent=4)
    return unique_variants


if __name__ == '__main__':
    get_variant()
    print('ok')
