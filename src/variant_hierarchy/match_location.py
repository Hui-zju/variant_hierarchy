from .get_accession import DEFAULT_UTA_DSN
from .get_n_interval import to_n_interval


def location_match(location_a, location_b, uta_dsn=DEFAULT_UTA_DSN):
    """
    Return the relationship between two variant locations:
      - "equal"
      - "contains"
      - "unrelated"
    """
    interval_a = to_n_interval(location_a, uta_dsn=uta_dsn)
    interval_b = to_n_interval(location_b, uta_dsn=uta_dsn)

    if interval_a and interval_b and None not in interval_a and None not in interval_b:
        set_a = set(interval_a)
        set_b = set(interval_b)

        if set_a == set_b:
            return "equal"

        if len(set_a) == 1 and set_a.issubset(set_b):
            return "contains"

        if len(set_a) == 1 and len(set_b) == 1:
            a = interval_a[0]
            b = interval_b[0]
            if a[0] == b[0] and a[1] <= b[1] and a[2] >= b[2]:
                return "contains"

    return "unrelated"


def test_location_match():
    location_a = {"gene": "EGFR"}
    location_b = {"gene": "EGFR"}
    assert location_match(location_a, location_b) == "equal"

    location_a = {"gene": "ALK", "gene2": "EGFR"}
    location_b = {"gene": "ALK", "gene2": "EGFR"}
    assert location_match(location_a, location_b) == "equal"

    location_a = {"gene": "ALK"}
    location_b = {"gene": "ALK", "gene2": "EGFR"}
    assert location_match(location_a, location_b) == "contains"

    location_a = {"gene": "EGFR"}
    location_b = {"gene": "EGFR", "exon": "3-5"}
    assert location_match(location_a, location_b) == "contains"

    location_a = {"gene": "EGFR", "exon": "3-5"}
    location_b = {"gene": "EGFR", "exon": "4"}
    assert location_match(location_a, location_b) == "contains"

    location_a = {"gene": "ALK"}
    location_b = {"gene": "EGFR"}
    assert location_match(location_a, location_b) == "unrelated"

    location_a = {"gene": "ALK", "gene2": "EGFR"}
    location_b = {"gene": "ALK"}
    assert location_match(location_a, location_b) == "unrelated"

    location_a = {"exon": "3"}
    location_b = {"gene": "EGFR"}
    assert location_match(location_a, location_b) == "unrelated"


if __name__ == "__main__":
    test_location_match()
