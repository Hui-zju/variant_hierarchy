import sys
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from variant_hierarchy.normalize.variant_parser import parse_variant


def norm_variant(gene, child_variant, parent_variant):

    parent = parse_variant(parent_variant, gene, variant_class=None)
    child = parse_variant(child_variant, gene, variant_class=None)
    if parent and child:
        # child['function'].update(parent['function'])
        variants = [child, parent]
        return variants
    if parent is None:
        print(f"Unable to parse variant from variantName (parent_variant={parent_variant}),gene={gene}")
    if child is None:
        print(f"Unable to parse variant from variantName (child_variant={child_variant}),gene={gene}")
    return []


if __name__ == '__main__':
    oncokb_function_variant_path = "../data/sources/Oncokb/oncokb_variant_function_annotation.csv"
    file_path = "/mnt/Data.Common/Chenh/PycharmProjects/InformationExtraction/dataset/graphkb_dataset/raw/Oncokb/oncokb_variant_function_annotation.csv"
    records = pd.read_csv(oncokb_function_variant_path).to_dict(orient='records')
    oncokb_function_variants = [parse_variant(record['variant'], record['gene'], variant_class=None) for record in records]

    cgi_function_variant_path = "../data/sources/CGI/catalog_of_validated_oncogenic_mutations_latest/catalog_of_validated_oncogenic_mutations.tsv"
    records = pd.read_csv(cgi_function_variant_path, sep='\t').to_dict(orient='records')
    cgi_function_variants = [parse_variant(record['protein'], record['gene'], variant_class=None) for record in records]

    PMC4232638_function_variant_path = "../data/sources/PMC4232638/13059_2014_484_MOESM2_ESM.csv"
    records = pd.read_csv(PMC4232638_function_variant_path).to_dict(orient='records')
    PMC4232638_function_variants = [parse_variant(record['Amino acid change'], record['gene'], variant_class=None) for record in records]





