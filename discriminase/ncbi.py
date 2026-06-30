"""NCBI Entrez access: list candidate genomes and fetch sequences (cached)."""
import os
import re

import certifi
from Bio import Entrez

os.environ.setdefault('SSL_CERT_FILE', certifi.where())


def setup_entrez(email: str) -> None:
    if not email:
        raise ValueError("NCBI requires an email. Set NCBI_EMAIL or pass --email.")
    Entrez.email = email


# --- cached raw-sequence fetches (used by the streaming builder) ------------
def _cache_path(cache_dir: str, key: str) -> str:
    safe = re.sub(r"[^A-Za-z0-9._-]+", "_", key)
    return os.path.join(cache_dir, f"{safe}.fasta")


def _read_fasta_text(text: str) -> str:
    return "".join(line.strip() for line in text.splitlines() if not line.startswith(">"))


def fetch_sequence_by_accession(accession: str, cache_dir: str = None):
    """Return (sequence, meta) for an NCBI nucleotide accession, caching the FASTA."""
    if cache_dir:
        os.makedirs(cache_dir, exist_ok=True)
        path = _cache_path(cache_dir, accession)
        if os.path.exists(path):
            with open(path, encoding="utf-8") as fh:
                return _read_fasta_text(fh.read()), {"accession": accession, "taxid": None}
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
    text = handle.read()
    handle.close()
    if cache_dir:
        with open(_cache_path(cache_dir, accession), "w", encoding="utf-8") as fh:
            fh.write(text)
    return _read_fasta_text(text), {"accession": accession, "taxid": None}


def fetch_sequence_by_taxid(taxid: str, cache_dir: str = None, name: str = None):
    """Return (sequence, meta) for a representative complete genome of ``taxid``.

    Picks the first complete genome for the taxon -- acceptable for a *commensal*
    panel, where any representative assembly is fine. Precise, reproducible picking
    matters for the *target* and is handled by the candidate-listing search (Phase 3).
    """
    if cache_dir:
        os.makedirs(cache_dir, exist_ok=True)
        path = _cache_path(cache_dir, f"txid{taxid}")
        if os.path.exists(path):
            with open(path, encoding="utf-8") as fh:
                return _read_fasta_text(fh.read()), {"accession": None, "taxid": str(taxid)}
    term = f"txid{taxid}[Organism:exp] AND complete genome[Title]"
    handle = Entrez.esearch(db="nucleotide", term=term, retmax=1)
    ids = Entrez.read(handle).get("IdList", [])
    handle.close()
    if not ids:
        raise ValueError(f"no complete genome found for taxid {taxid}")
    handle = Entrez.efetch(db="nucleotide", id=ids[0], rettype="fasta", retmode="text")
    text = handle.read()
    handle.close()
    accession = text[1:].split()[0] if text.startswith(">") else None
    if cache_dir:
        with open(_cache_path(cache_dir, f"txid{taxid}"), "w", encoding="utf-8") as fh:
            fh.write(text)
    return _read_fasta_text(text), {"accession": accession, "taxid": str(taxid)}


def search_target_candidates(query: str, retmax: int = 25, refseq_only: bool = True):
    """List candidate genome assemblies for a name/term -- never auto-pick one.

    A species name maps to many assemblies; the reproducible key is the accession.
    This returns enough to choose from: accession, full title, organism, length,
    taxid. The caller (CLI prompt or web table) picks; nothing is selected silently.
    """
    parts = [f"{query}[Organism] OR {query}[Title]", "complete genome[Title]"]
    if refseq_only:
        parts.append("srcdb_refseq[PROP]")
    term = " AND ".join(f"({p})" for p in parts)
    handle = Entrez.esearch(db="nucleotide", term=term, retmax=retmax)
    ids = Entrez.read(handle).get("IdList", [])
    handle.close()
    if not ids:
        return []
    handle = Entrez.esummary(db="nucleotide", id=",".join(ids))
    summaries = Entrez.read(handle)
    handle.close()
    out = []
    for s in summaries:
        out.append(
            {
                "accession": s.get("AccessionVersion") or s.get("Caption"),
                "title": str(s.get("Title", "")),
                "organism": str(s.get("Organism", "")) or _organism_from_title(str(s.get("Title", ""))),
                "length_bp": int(s.get("Length") or s.get("Slen") or 0),
                "taxid": str(int(s.get("TaxId") or 0)),   # TaxId is an Entrez IntegerElement
            }
        )
    out.sort(key=lambda r: r["length_bp"], reverse=True)
    return out


def _organism_from_title(title: str) -> str:
    # "Salmonella enterica subsp. ... chromosome, complete genome" -> leading binomial
    words = title.split()
    return " ".join(words[:2]) if len(words) >= 2 else title
