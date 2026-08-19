import json
import sys
import pandas as pd
from pathlib import Path
from copy import deepcopy

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from variant_norm_civic import norm_variant as norm_civic
from variant_norm_cgi import norm_variant as norm_cgi
from variant_norm_oncokb import norm_variant as norm_oncokb
from variant_norm_cosmic import norm_variant as norm_cosmic
from variant_norm_function import norm_variant as norm_function

CONFIG = [
    {
        "source": "CIVIC",
        "file_path": "../data/sources/CIVIC/01-Nov-2024-VariantSummaries.tsv",
        "norm_function": norm_civic,
        "function_parameters": ["Variable.variant", "Variable.feature_name", "Variable.variant_types"],
        "is_dedup": True,
        "save_path": "../data/variants/civic_variants.json",
    },
    {
        "source": "OncoKB",
        "file_path": "../data/sources/Oncokb/oncokb_biomarker_drug_associations.tsv",
        "norm_function": norm_oncokb,
        "function_parameters": ["Variable.Alterations", "Variable.Gene"],
        "is_dedup": True,
        "save_path": "../data/variants/oncokb_variants.json",
    },
    {
        "source": "COSMIC",
        "file_path": "../data/sources/Cosmic/Actionability_AllData_v13_GRCh37.tsv",
        "norm_function": norm_cosmic,
        "function_parameters": ["Variable.MUTATION_REMARK"],
        "is_dedup": True,
        "save_path": "../data/variants/cosmic_variants.json",
    },
    {
        "source": "CGI",
        "file_path": "../data/sources/CGI/cancergenomeinterpreter.org_data_biomarkers_cgi_biomarkers_latest.tsv",
        "norm_function": norm_cgi,
        "function_parameters": ["Variable.Biomarker", "Variable.Alteration type"],
        "is_dedup": True,
        "save_path": "../data/variants/cgi_variants.json",
    },
    {
        "source": "OncoKB",
        "file_path": "../data/sources/Oncokb/oncokb_variant_function_annotation.csv",
        "norm_function": norm_function,
        "function_parameters": ["Variable.gene", "Variable.variant", "Variable.knownEffect"],
        "is_dedup": True,
        "save_path": "../data/variants/oncokb_function_variants.json",
    },
    {
        "source": "CGI",
        "file_path": "../data/sources/CGI/catalog_of_validated_oncogenic_mutations_latest/catalog_of_validated_oncogenic_mutations.tsv",
        "norm_function": norm_function,
        "function_parameters": ["Variable.gene", "Variable.protein", "Fixed.oncogenic"],
        "is_dedup": True,
        "save_path": "../data/variants/cgi_function_variants.json",
    },
    {
        "source": "PMC4232638",
        "file_path": "../data/sources/PMC4232638/13059_2014_484_MOESM2_ESM.csv",
        "norm_function": norm_function,
        "function_parameters": ["Variable.Gene", "Variable.Amino acid change", "Variable.Type"],
        "is_dedup": True,
        "save_path": "../data/variants/PMC4232638_function_variants.json",
    }
]


def resolve_config_path(path):
    path = Path(path)
    if path.is_absolute():
        return path
    return (Path(__file__).resolve().parent / path).resolve()


def convert_set_to_list(obj):
    if isinstance(obj, set):
        return list(obj)
    elif isinstance(obj, dict):
        return {key: convert_set_to_list(value) for key, value in obj.items()}
    elif isinstance(obj, list):
        return [convert_set_to_list(item) for item in obj]
    else:
        return obj

def get_variants(source, file_path, norm_function, function_parameters, is_dedup=True, save_path=None):
    file_path = resolve_config_path(file_path)
    if file_path.suffix == ".tsv":
        records = pd.read_csv(file_path, sep='\t').to_dict(orient='records')
    elif file_path.suffix == ".csv":
        records = pd.read_csv(file_path).to_dict(orient='records')
    else:
        raise ValueError(f"unsupported file extension: {file_path}")

    variants = []
    for record in records:
        args = []
        raw = {}
        for key in function_parameters:
            if key.startswith("Variable."):
                new_key = key.removeprefix("Variable.")
                v = record.get(new_key)
            elif key.startswith("Fixed."):
                new_key = key.removeprefix("Fixed.")
                v = new_key
            else:
                raise ValueError("unknown key {}".format(key))
            args.append(str(v))
            raw[new_key] = v

        items = norm_function(*args)
        for item in items:
            variants.append({
                **item,
                "source": [{source: raw}],
            })

    if is_dedup:
        variants = list({json.dumps(v, sort_keys=True): v for v in variants}.values())

    if save_path:
        save_path = resolve_config_path(save_path)
        save_path.parent.mkdir(parents=True, exist_ok=True)
        with open(save_path, "w", encoding="utf-8") as f:  # ../data/variants/
            json.dump(variants, f, ensure_ascii=False, indent=4)

    return variants


def run_from_config(config_list):
    results = []
    for c in config_list:
        print(f"--------------------------START: {c['file_path']}--------------------------------")
        variants = get_variants(c["source"],
            c["file_path"],
            c["norm_function"],
            c["function_parameters"],
            c["is_dedup"],
            c["save_path"]
        )
        results += variants

    merged_dict = {}
    for item in results:
        k = item['display_name']
        if k not in merged_dict:
            merged_dict[k] = deepcopy(item)
        else:
            m = merged_dict[k]
            assert {k:v for k,v in m.items() if k not in ('source','function')} == {k:v for k,v in item.items() if k not in ('source','function')}
            m['source'].extend(item['source'])
            m['function'].update(item['function'])

    all_variants = list(merged_dict.values())
    all_variants = [convert_set_to_list(item) for item in all_variants]

    return all_variants


def get_variants_with_streamlined_source(variants):
    variants_with_streamlined_source = []
    for var in variants:
        source = var['source']
        streamlined_source = sorted({next(iter(item)) for item in source if item})
        new_var = var.copy()
        new_var['source'] = streamlined_source

        variants_with_streamlined_source.append(new_var)
    return variants_with_streamlined_source


def save_variants(variants, save_path):
    save_path = resolve_config_path(save_path)
    save_path.parent.mkdir(parents=True, exist_ok=True)
    with open(save_path, "w", encoding="utf-8") as f:
        json.dump(variants, f, ensure_ascii=False, indent=4)


if __name__ == '__main__':
    save_path = "../data/variants/variants.json"
    variants = run_from_config(CONFIG)
    variants = get_variants_with_streamlined_source(variants)
    save_variants(variants, save_path)
    print('ok')



