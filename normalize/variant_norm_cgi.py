import re
import json
import sys
import pandas as pd
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from variant_hierarchy.normalize.variant_parser import parse_variant
from source_mappings.loader import load_yaml_mapping
CGI_MAPPING = load_yaml_mapping("cgi_biomarker_aliases.yaml")
CGI_VARIANT_TYPE_TO_CLASS = load_yaml_mapping("cgi_alteration_type_classes.yaml")


def split_cgi_variant_types(variant_class):
    value = str(variant_class or "").strip().upper()
    if not value:
        return set()
    return {item.strip() for item in value.split(";") if item.strip()}


def infer_variant_class_from_cgi(variant_class):
    variant_types = split_cgi_variant_types(variant_class)
    if len(variant_types) > 1:
        return None

    classes = {
        CGI_VARIANT_TYPE_TO_CLASS[item]
        for item in variant_types
        if item in CGI_VARIANT_TYPE_TO_CLASS
    }
    return classes.pop() if len(classes) == 1 else None


def norm_variant(biomarker, variant_class):
    biomarker = CGI_MAPPING.get(biomarker.strip(), biomarker.strip())
    normalized_variant_class = infer_variant_class_from_cgi(variant_class)


    # ----------------------------------------split variant------------------------------------------
    # 'AR (F877L) + AR (T878A)'
    if len(re.split(r'\s*\+\s*(?!$)', biomarker)) > 1:
        result = []
        for n in re.split(r'\s*\+\s*', biomarker):
            result.extend(norm_variant(n.strip(), variant_class))
        return result

    # ABL1 (F359V,F359C,F359I,Y253H,E255K,E255V)  JAK1 (S646F;R683)
    match = re.match(r"^(\w+) \(([A-Z0-9*.,;\s]+)\)$", biomarker)
    if match:
        gene = match.group(1)
        mutations = re.split(r'[;,]', match.group(2)) # match.group(2).split(",")
        variants = [f"{gene} {mutation.strip()}" for mutation in mutations]
        result = []
        for n in variants:
            result.extend(norm_variant(n.strip(), variant_class))
        return result

    # KIT mutation in exon 9,11,13,14 or 17
    match = re.match(r'(\w+)\s+mutation\s+in\s+exon\s*(\d+(?:,\d+)*)\s*(?:or)?\s*(\d+)*', biomarker)
    if match:
        result = []
        gene = match.group(1)
        exons = match.group(2)
        exon_numbers = exons.split(',')
        variants = [f"{gene} exon {exon.strip()}" for exon in exon_numbers]
        for variant in variants:
            result.extend(norm_variant(variant, variant_class))
        return result

    # PDGFRA (552-596,631-668,814-854)
    match = re.match(r'^(\w+)\s\(([\d,-]+)\)$', biomarker)
    if match and ',' in biomarker:
        result = []
        gene = match.group(1)
        ranges = match.group(2)
        range_list = ranges.split(',')
        variants = [f"{gene} ({r})" for r in range_list]
        for variant in variants:
            result.extend(norm_variant(variant, variant_class))
        return result


    # CSF3R frameshift variant (D771),frameshift variant (S783)
    match = re.match(r"^([\w\d]+)\s(.+)$", biomarker)
    if match and ',' in biomarker:
        result = []
        gene_name = match.group(1)
        variant_name = match.group(2)

        variant_parts = re.split(r"\s*,\s*", variant_name)
        variants = [f"{gene_name} {variant}" for variant in variant_parts]
        for variant in variants:
            result.extend(norm_variant(variant, variant_class))
        return result

    # ----------------------------------------norm variant------------------------------------------


    parsed_variant = parse_variant(biomarker, variant_class=normalized_variant_class)
    if parsed_variant:
        return [parsed_variant]
    print(f"Unable to parse variant from variantName (variantName={biomarker})")
    return []
    # raise Exception(f"Unable to parse variant from variantName (variantName={biomarker})")



def preprocess_rows(rows):
    new_rows = []
    for row in rows:
        biomarker = row.get('Biomarker', '')
        variant_class = row.get('Alteration type', '')
        proteins = norm_variant(biomarker, variant_class)
        for protein in proteins:
            new_row = row.copy()
            new_row['protein'] = protein
            new_row['protein']['Alteration'] = row['Alteration']
            new_row['protein']['Biomarker'] = row['Biomarker']
            new_row['protein']['Alteration type'] = row['Alteration type']
            new_rows.append(new_row)
    return new_rows


def get_variant():
    file_path = "../data/sources/CGI/cancergenomeinterpreter.org_data_biomarkers_cgi_biomarkers_latest.tsv"
    df = pd.read_csv(file_path, sep='\t')
    records = df.to_dict(orient='records')

    variants = []
    for record in records:
        items = norm_variant(record['Biomarker'], record['Alteration type'])
        for item in items:
            variants.append({**item, **{"source": ["CGI"]},
                             **{"raw": {"alteration": record['Biomarker'], "alteration_type": record['Alteration type']}}})


    unique_variants = list({json.dumps(v, sort_keys=True): v for v in variants}.values())
    # gene_variants = [variant for variant in unique_variants if variant['variantClass'] == "GeneVariant"]
    # unresolved_variants = [variant for variant in unique_variants if variant['variantClass'] == "UnresolvedVariant"]
    with open("cgi_variants.json", "w", encoding="utf-8") as f: # ../data/variants/
        json.dump(unique_variants, f, ensure_ascii=False, indent=4)
    return unique_variants


if __name__ == '__main__':
    # variants = get_variant()
    # variants = [item for variant in variants for item in variant]
    # dic = {'aa_change': '61', 'displayName': 'KRAS 61', 'reference1': 'KRAS'}
    # if dic in variants:
    #     print('error')

    # res = norm_variant('ERBB2 inframe insertion (P780GSP),inframe insertion (781GSP),inframe insertion (A775YVMA),inframe insertion (G776YVMA)', '')

    # typ = df['Alteration type'].unique()
    # rows = preprocess_rows(records)
    # variants = [row['protein'] for row in rows]
    # variants = [item for variant in variants for item in variant]
    # test(variants)

    res = parse_variant("EGFR oncogenic")
    if res:
        print(res)

    print("ok")
