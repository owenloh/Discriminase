"""Screen a target's guides against the commensal index.

A target guide is kept (commensal-sparing) when the index finds no commensal cut-site
within the seed-anchored distance model -- i.e. :meth:`GuideIndex.query` returns
``None``. GC filtering is applied to the target guides first (it only makes sense for
the guides we might actually use, not for the commensals we are sparing).
"""
import sys
from typing import Callable, List, Optional

import numpy as np

from .config import GuideFinderConfig
from .genome import extract_target_guides, gc_fraction
from .index import GuideIndex
from .pack import unpack_guide


def _progress(done: int, total: int, kept: int) -> None:
    if total == 0:
        return
    bar = "#" * int(24 * done / total)
    sys.stdout.write(f"\r  screening [{bar:<24}] {done}/{total}  ({kept} kept)")
    sys.stdout.flush()


def find_sparing_guides(
    target_seq: str,
    index: GuideIndex,
    config: GuideFinderConfig,
    progress: bool = True,
    on_progress: Optional[Callable[[int, int, int], None]] = None,
) -> List[dict]:
    """Return commensal-sparing guides for ``target_seq`` as a list of row dicts."""
    L = config.guide_length
    packed, starts, strands = extract_target_guides(
        target_seq, L, config.pam, config.pam_to_guide_gap, config.pam_side
    )
    flip = config.pam_side == "3prime"   # stored PAM-proximal-first -> show 5'->3'

    # GC filter (target guides only).
    if packed.size:
        frac = gc_fraction(packed, L)
        keep = (frac >= config.min_gc) & (frac <= config.max_gc)
        packed, starts, strands, frac = packed[keep], starts[keep], strands[keep], frac[keep]
    else:
        frac = np.empty(0)

    total = int(packed.size)
    d, s = config.total_max_mismatch, config.seed_max_mismatch
    cap = config.max_guides or total
    organisms = index.organisms

    rows: List[dict] = []
    report = on_progress or (_progress if progress else None)
    for i in range(total):
        if report and i % 500 == 0:
            report(i, total, len(rows))
        g = int(packed[i])
        if index.query(g, d, s) is not None:
            continue                       # collides with a commensal -> drop
        seq = unpack_guide(g, L)
        rows.append(
            {
                "position": int(starts[i]),
                "strand": "-" if strands[i] else "+",
                "guide_sequence": seq[::-1] if flip else seq,
                "gc": round(float(frac[i]), 3),
            }
        )
        if len(rows) >= cap:
            break

    if report:
        report(total, total, len(rows))
        if progress and not on_progress:
            sys.stdout.write("\n")
    return rows
