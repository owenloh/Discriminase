"""Central configuration for Discriminase.

Every tunable lives here so the CLI, the database builder, the search engine and the
web UI all read from one place. Construct a ``GuideFinderConfig()`` and override
fields as needed.
"""
import os
from typing import List, Optional

from pydantic import BaseModel


def _default_workers() -> int:
    """Use the machine's cores, leaving a little headroom; never block the box."""
    cpus = os.cpu_count() or 1
    return max(1, cpus - 1)


# Convenience starting points; users can still override any field. A preset only
# sets the nuclease geometry (PAM, side, typical guide/seed lengths).
NUCLEASE_PRESETS = {
    "cas12a": {"pam": "TTTV", "pam_side": "5prime",
               "guide_length": 23, "seed_len": 10},
    "spcas9": {"pam": "NGG", "pam_side": "3prime",
               "guide_length": 20, "seed_len": 10},
}


class GuideFinderConfig(BaseModel):
    # --- guide geometry -----------------------------------------------------
    guide_length: int = 23          # spacer length in nt (must be <= 31 for uint64)

    # --- PAM / nuclease -----------------------------------------------------
    # PAM accepts IUPAC codes (N, R, Y, V, ...). `pam_side` is "5prime" (Cas12a:
    # PAM before the spacer) or "3prime" (SpCas9 NGG: PAM after the spacer); the
    # seed is always the PAM-proximal end, whichever side that is. The protospacer
    # is adjacent to the PAM -- to leave room between them, pad the PAM with N
    # (e.g. "TTTN"). See NUCLEASE_PRESETS.
    pam: str = "TTT"
    pam_side: str = "5prime"

    # --- distance model (seed-anchored; see docs/ARCHITECTURE.md) ------------
    # A target guide collides with a commensal cut-site when their PAM-proximal
    # `seed_len` bases match within `seed_max_mismatch` AND the full guides match
    # within `total_max_mismatch`. Defaults are conservative-but-reasonable and
    # tunable per nuclease -- they are a model, not dogma.
    seed_len: int = 10
    seed_max_mismatch: int = 1
    total_max_mismatch: int = 4

    # --- target guide filters (applied to the TARGET only, not commensals) --
    min_gc: float = 0.4
    max_gc: float = 0.6
    max_guides: int = 1000           # stop after this many surviving guides (0 = all)

    # --- NCBI Entrez --------------------------------------------------------
    # NCBI asks callers to identify themselves: set NCBI_EMAIL or pass --email.
    email: str = os.environ.get("NCBI_EMAIL", "")

    # --- compute ------------------------------------------------------------
    # 0 -> adapt to the machine (os.cpu_count); set a number to pin it. This is
    # read through `n_workers` so changing hardware needs no config edit.
    processes: int = 0

    # --- database location --------------------------------------------------
    # Files: "{db_prefix}_len{guide_length}.idx.npy / .org.npy / .meta.json"
    db_dir: str = "database"
    db_prefix: str = "protected_guides"

    # --- genome cache (downloaded FASTAs, reused across builds) -------------
    cache_dir: str = os.path.join("database", "genomes")

    @property
    def n_workers(self) -> int:
        return self.processes if self.processes and self.processes > 0 else _default_workers()

    def index_prefix(self) -> str:
        return os.path.join(self.db_dir, f"{self.db_prefix}_len{self.guide_length}")
