"""Algorithm selection that travels with channel weights into reduction workers."""

from dataclasses import dataclass


def normalize_demosaic_algorithm(value):
    algorithm = str(value or 'bilinear').strip().lower() or 'bilinear'
    if algorithm == 'ddfapd':
        algorithm = 'menon2007'
    if algorithm not in {'bilinear', 'malvar2004', 'menon2007'}:
        raise ValueError(
            f'Invalid Demosaic Algorithm {value!r}; choose bilinear, malvar2004, or menon2007 (ddfapd).'
        )
    return algorithm


@dataclass(frozen=True)
class DemosaicMix:
    """Picklable non-default algorithm and normalized output-channel weights."""

    weights: object
    algorithm: str


def reconstruct_bayer(image, pattern, algorithm):
    # Lazy lookup retains compatibility with installations/tests using only bilinear.
    import colour_demosaicing

    name = {'bilinear': 'demosaicing_CFA_Bayer_bilinear',
            'malvar2004': 'demosaicing_CFA_Bayer_Malvar2004',
            'menon2007': 'demosaicing_CFA_Bayer_Menon2007'}[normalize_demosaic_algorithm(algorithm)]
    return getattr(colour_demosaicing, name)(image, pattern)
