"""Central configuration for Discriminase.

Every tunable lives here so the CLI, the database builder and the search engine
all read from one place. Construct a ``GuideFinderConfig()`` and override fields
as needed (the CLI does exactly this from its flags).
"""
import os
from typing import List, Optional

from pydantic import BaseModel


class GuideFinderConfig(BaseModel):
    # --- guide geometry -----------------------------------------------------
    guide_length: int = 23          # spacer length in nucleotides
    seed_len: int = 4               # reserved for seed-region scoring

    # --- PAM (Cas12a / Cpf1 style: a 5' T-rich PAM upstream of the spacer) --
    # NOTE: this is a Cas12a-style 5' PAM, *not* SpCas9 (NGG, 3'). Matched
    # exactly as written, so "TTT" also accepts a following T. Change `pam`
    # here if you target a different nuclease.
    pam: str = "TTT"
    pam_to_guide_gap: int = 1       # bases between PAM and spacer start

    # --- NCBI Entrez --------------------------------------------------------
    # NCBI asks every caller to identify themselves. Set the NCBI_EMAIL env var
    # or pass --email on the CLI; this default is only a fallback.
    email: str = os.environ.get("NCBI_EMAIL", "owenloh0607@gmail.com")

    # --- encoding -----------------------------------------------------------
    all_N_into_A: bool = True       # map ambiguous N bases to A (see encoding.py)

    # --- guide filters ------------------------------------------------------
    min_gc: float = 0.4
    max_gc: float = 0.6
    similarity_threshold: float = 0.70  # max_mismatches = int((1 - this) * L)
    if_sig_cutoff: bool = True          # use the trie seed prefilter
    sig_cutoff: int = 13                # exact-match seed length for the trie

    # --- search controls ----------------------------------------------------
    subsection: List[Optional[int]] = [0, None]  # slice of target guides to scan
    max_guides: int = 1000                       # stop after this many hits
    processes: int = 4                           # PAM scan worker processes

    # --- database location --------------------------------------------------
    # Files are named "{db_prefix}_len{guide_length}.trie/.bktree" under db_dir.
    db_dir: str = "database"
    db_prefix: str = "protected_guides"

    @property
    def max_mismatches(self) -> int:
        return int((1 - self.similarity_threshold) * self.guide_length)

    def trie_path(self) -> str:
        return os.path.join(self.db_dir, f"{self.db_prefix}_len{self.guide_length}.trie")

    def bktree_path(self) -> str:
        return os.path.join(self.db_dir, f"{self.db_prefix}_len{self.guide_length}.bktree")
