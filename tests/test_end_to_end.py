"""Build -> find on controlled synthetic genomes (no network, no external files).

Pins the whole pipeline AND the safety contract: a target guide that is also a
commensal cut-site must be dropped; a genuinely unique guide must survive.
"""
import numpy as np

from discriminase.config import GuideFinderConfig
from discriminase.database import build_index
from discriminase.search import find_sparing_guides

# 23-mers, ~50% GC so they pass the GC filter. SHARED is planted in both the
# commensal and the target; UNIQUE only in the target and far from everything.
SHARED = "ACGTACGTACGTACGTACGTACG"
UNIQUE = "TGCATGCATGCATGCATGCATGC"
PAD = "GA" * 12                       # no "TTT", so it spawns no forward PAM sites


def _cassette(spacer):                # 5' TTT PAM, then the adjacent spacer
    return "TTT" + spacer


def _cfg(tmp_path):
    cfg = GuideFinderConfig()
    cfg.db_dir = str(tmp_path)
    cfg.db_prefix = "panel"
    cfg.email = ""
    cfg.max_guides = 0                 # keep all survivors
    return cfg


def test_build_then_find(tmp_path):
    cfg = _cfg(tmp_path)
    commensal = PAD + _cassette(SHARED) + PAD
    target = PAD + _cassette(SHARED) + PAD + _cassette(UNIQUE) + PAD

    index = build_index(
        [{"name": "synthetic-commensal", "seq": commensal}],
        cfg,
        get_sequence=lambda src: (src["seq"], {}),
        progress=lambda *_: None,
    )
    assert len(index) > 0

    rows = find_sparing_guides(target, index, cfg, progress=False)
    kept = {r["guide_sequence"] for r in rows}

    assert UNIQUE in kept, "a genuinely unique guide must survive"
    assert SHARED not in kept, "a guide identical to a commensal must be dropped (unsafe!)"


def test_self_panel_keeps_nothing(tmp_path):
    """Screening a genome against an index built from itself yields no guides."""
    cfg = _cfg(tmp_path)
    rng = np.random.default_rng(3)
    genome = "".join(rng.choice(list("ACGT"), 5000))

    index = build_index(
        [{"name": "self", "seq": genome}],
        cfg,
        get_sequence=lambda src: (src["seq"], {}),
        progress=lambda *_: None,
    )
    rows = find_sparing_guides(genome, index, cfg, progress=False)
    assert rows == []
