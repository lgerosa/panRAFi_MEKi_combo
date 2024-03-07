#Return the integer species number for a given species in the model

from pysb.core import as_complex_pattern

def get_species_index(model, pattern):
    pattern = as_complex_pattern(pattern)
    matches = [
        i for i, s in enumerate(model.species)
        if s.is_equivalent_to(pattern)
    ]
    n = len(matches)
    assert n == 1, f"Expected exactly one match, got {n}"
    return matches[0]

