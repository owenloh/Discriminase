"""PAM anchoring and guide (spacer) extraction.

A guide is a fixed-length window adjacent to a PAM. We scan every position of a
genome on both strands, record the spacer start wherever a PAM matches, then
keep spacers whose GC content is in range. Genome scanning parallelises across
genomes; for a single genome it runs serially.
"""
from functools import partial
from typing import List, Tuple
from concurrent.futures import ProcessPoolExecutor, as_completed

from bitarray import bitarray

from .config import GuideFinderConfig
from .encoding import TwoBitDNA


def process_genome(genome: TwoBitDNA, pam_bits: bitarray, pam_length: int, gap: int):
    """Return [(genome, spacer_start, is_rev), ...] for every PAM hit in `genome`."""
    guides: List[Tuple[TwoBitDNA, int, bool]] = []
    L = len(genome)
    buffer_length = pam_length + gap

    for i in range(L):
        forward_start = i + buffer_length
        backward_start = i - buffer_length

        if forward_start < L and genome.get_encoded_slice(i, pam_length, False) == pam_bits:
            guides.append((genome, forward_start, False))

        if backward_start >= buffer_length and genome.get_encoded_slice(i, pam_length, True) == pam_bits:
            guides.append((genome, backward_start, True))

    return guides


def build_pam_valid_guides(genomes: List[TwoBitDNA], config: GuideFinderConfig):
    """All PAM-anchored spacer positions across `genomes`."""
    pam = TwoBitDNA("PAM", config.pam)
    worker = partial(process_genome, pam_bits=pam.bits, pam_length=pam.length, gap=config.pam_to_guide_gap)

    guides: List[Tuple[TwoBitDNA, int, bool]] = []
    if config.processes <= 1 or len(genomes) <= 1:
        for g in genomes:
            guides.extend(worker(g))
        return guides

    with ProcessPoolExecutor(max_workers=config.processes) as executor:
        futures = [executor.submit(worker, g) for g in genomes]
        for future in as_completed(futures):
            guides.extend(future.result())
    return guides


def gc_count(bits: bitarray) -> int:
    """Count G/C bases. C=01, G=10 -> the two bits differ; A=00, T=11 -> equal."""
    gc = 0
    for i in range(0, len(bits), 2):
        if bits[i] != bits[i + 1]:
            gc += 1
    return gc


def build_valid_target_guides(genome: TwoBitDNA, config: GuideFinderConfig):
    """Extract GC-filtered candidate guides from a single target genome.

    Returns [(spacer_start, is_rev, bits), ...].
    """
    valid: List[Tuple[int, bool, bitarray]] = []
    length = config.guide_length

    for _, index, is_rev in build_pam_valid_guides([genome], config):
        bits = genome.get_encoded_slice(index, length, is_rev)
        if len(bits) != 2 * length:
            continue  # window ran off the end of the genome
        if config.min_gc <= gc_count(bits) / length <= config.max_gc:
            valid.append((index, is_rev, bits))

    return valid
