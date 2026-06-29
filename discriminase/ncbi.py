"""NCBI Entrez access: resolve organism names to taxids and fetch genomes."""
import os

import certifi
from Bio import Entrez, SeqIO

from .encoding import TwoBitDNA

os.environ.setdefault('SSL_CERT_FILE', certifi.where())


def setup_entrez(email: str) -> None:
    if not email:
        raise ValueError("NCBI requires an email. Set NCBI_EMAIL or pass --email.")
    Entrez.email = email


def resolve_taxid(name: str):
    """Resolve an organism name to (taxid, scientific_name), or None if not found."""
    handle = Entrez.esearch(db="taxonomy", term=name, retmax=1)
    record = Entrez.read(handle)
    handle.close()
    ids = record.get("IdList", [])
    if not ids:
        return None
    taxid = ids[0]
    handle = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
    data = Entrez.read(handle)
    handle.close()
    sci_name = data[0]["ScientificName"] if data else name
    return taxid, sci_name


def _fetch_first_genome(term: str, label: str):
    handle = Entrez.esearch(db="nucleotide", term=term, retmax=5)
    record = Entrez.read(handle)
    handle.close()
    ids = record.get("IdList", [])
    if not ids:
        return None
    seq_record = SeqIO.read(
        Entrez.efetch(db="nucleotide", id=ids[0], rettype="fasta", retmode="text"),
        "fasta",
    )
    return TwoBitDNA(label, str(seq_record.seq))


def fetch_genome_by_name(name: str) -> TwoBitDNA:
    """Fetch a complete RefSeq genome by organism name."""
    term = f'{name}[Organism] AND "complete genome"[Title] AND srcdb_refseq[PROP]'
    genome = _fetch_first_genome(term, name)
    if genome is None:
        raise ValueError(f"No complete RefSeq genome found for: {name}")
    return genome


def fetch_genome_by_taxid(taxid: str, name: str = None) -> TwoBitDNA:
    """Fetch a complete genome by NCBI taxonomy id."""
    label = name or f"txid{taxid}"
    term = f"txid{taxid}[Organism:exp] AND complete genome[Title]"
    genome = _fetch_first_genome(term, label)
    if genome is None:
        raise ValueError(f"No complete genome found for taxid {taxid}")
    return genome
