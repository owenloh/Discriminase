"""Readers for input files: commensal panels (CSV) and target sequences."""
import csv
import os
from typing import List, Tuple


def read_taxid_csv(file_path: str) -> List[Tuple[str, str]]:
    """Read a panel CSV with headers `taxid` and `organism_strain`.

    Returns [(organism_strain, taxid), ...].
    """
    organisms: List[Tuple[str, str]] = []
    with open(file_path, 'r', encoding='utf-8', newline='') as fh:
        reader = csv.DictReader(fh)
        required = {'taxid', 'organism_strain'}
        if not reader.fieldnames or not required.issubset(set(reader.fieldnames)):
            raise ValueError(
                f"CSV must have headers {sorted(required)}; found {reader.fieldnames}"
            )
        for row in reader:
            taxid = (row.get('taxid') or '').strip()
            strain = (row.get('organism_strain') or '').strip()
            if taxid and strain:
                organisms.append((strain, taxid))
    return organisms


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
