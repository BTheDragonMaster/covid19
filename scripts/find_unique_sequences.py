#!/usr/bin/nv python

import os
from sys import argv

from parsers import parse_fasta
from writers import write_fasta, write_code_to_accession

dirname = os.path.dirname(__file__)
CODE_DIR = os.path.join(dirname, '../data/code_to_accesion.fasta')
UNIQUE_SEQ_DIR = os.path.join(dirname, '../data/unique_sequences.fasta')


def parse_fasta_id(fasta_id):
    fasta_id = fasta_id.split('|')[0]
    fasta_id = fasta_id.strip()
    return fasta_id

def reverse_fasta_dict(fasta_dict):
    sequence_to_id = {}
    for fasta_id, sequence in fasta_dict.items():
        fasta_id = parse_fasta_id(fasta_id)
        if not sequence in sequence_to_id:
            sequence_to_id[sequence] = []

        sequence_to_id[sequence].append(fasta_id)

    return sequence_to_id

def assign_code(sequence_to_id):
    code_to_sequence = {}
    code_to_accession = {}
    for i, (sequence, accessions) in enumerate(sequence_to_id.items()):
        code = 'seq_%.4d' % i
        code_to_sequence[code] = sequence
        code_to_accession[code] = accessions

    return code_to_sequence, code_to_accession

if __name__ == "__main__":
    fasta = argv[1]

    id_to_sequence = parse_fasta(fasta)
    sequence_to_id = reverse_fasta_dict(id_to_sequence)
    code_to_sequence, code_to_accession = assign_code(sequence_to_id)
    write_fasta(code_to_sequence, UNIQUE_SEQ_DIR)
    write_code_to_accession(code_to_accession, CODE_DIR)
    





