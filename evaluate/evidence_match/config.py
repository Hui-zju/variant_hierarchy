from __future__ import annotations

from evaluate.common import PROJECT_ROOT


PATIENT_DIR = PROJECT_ROOT / "data" / "patient"
GENIE_DIR = PATIENT_DIR / "genie_mskcc_processed"
MSK_DIR = PATIENT_DIR / "msk_impact_2017"
VARIANTS_PATH = PROJECT_ROOT / "data" / "variants" / "variants.json"
HIERARCHY_PATH = PROJECT_ROOT / "data" / "hierarchy" / "hierarchy.json"
RESULT_DIR = PROJECT_ROOT / "evaluate" / "evidence_match" / "results"
INTERMEDIATE_DIR = PROJECT_ROOT / "evaluate" / "evidence_match" / "intermediate"
GENIE_NORMALIZED_PATH = INTERMEDIATE_DIR / "genie_normalized_variants.tsv"
MSK_NORMALIZED_PATH = INTERMEDIATE_DIR / "msk_impact_2017_normalized_variants.tsv"

MSK_MUTATION_API_PATH = MSK_DIR / "data_mutations_api.tsv"
MSK_MUTATION_API_JSON_PATH = MSK_DIR / "data_mutations_api_full.json"
CBIOPORTAL_MUTATION_URL = (
    "https://www.cbioportal.org/api/molecular-profiles/"
    "msk_impact_2017_mutations/mutations/fetch?projection=DETAILED"
)

SOURCES = ["OncoKB", "CIVIC", "COSMIC", "CGI"]
