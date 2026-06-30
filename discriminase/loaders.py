"""Readers for input files: commensal panels (CSV) and target sequences."""
import csv
import os
from typing import List, Tuple


def read_panel_csv(file_path: str) -> List[Tuple[str, str, str]]:
    """Read a panel CSV. Needs `organism_strain` and at least one of `taxid`/`accession`.

    Accession is the more reliable key (it pins one exact complete genome); taxid works
    too. Returns [(organism_strain, taxid, accession), ...] (missing fields as "").
    """
    rows: List[Tuple[str, str, str]] = []
    with open(file_path, 'r', encoding='utf-8', newline='') as fh:
        reader = csv.DictReader(fh)
        fields = set(reader.fieldnames or [])
        if 'organism_strain' not in fields or not (fields & {'taxid', 'accession'}):
            raise ValueError(
                "CSV must have 'organism_strain' and at least one of 'taxid'/'accession'; "
                f"found {reader.fieldnames}"
            )
        for row in reader:
            name = (row.get('organism_strain') or '').strip()
            taxid = (row.get('taxid') or '').strip()
            accession = (row.get('accession') or '').strip()
            if name and (taxid or accession):
                rows.append((name, taxid, accession))
    return rows


def read_sequence_file(file_path: str) -> str:
    """Read a DNA sequence from .txt/.fa/.fasta/.fna (FASTA headers stripped)."""
    ext = os.path.splitext(file_path)[1].lower()
    with open(file_path, 'r', encoding='utf-8') as fh:
        lines = fh.readlines()
    if ext in ('.fa', '.fasta', '.fna'):
        seq = ''.join(line.strip() for line in lines if not line.startswith('>'))
    elif ext == '.txt':
        seq = ''.join(lines)
    else:
        raise ValueError(f"Unsupported extension {ext}; use .txt/.fa/.fasta/.fna")
    return seq.strip().replace('\n', '').replace('\r', '').replace(' ', '')
