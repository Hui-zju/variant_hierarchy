from pathlib import Path
import yaml


def load_yaml_mapping(filename: str) -> dict:
    mapping_path = (
        Path(__file__).resolve().parent
        / filename
    )

    with open(mapping_path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)