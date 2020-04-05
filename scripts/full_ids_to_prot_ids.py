#!/usr/bin/env python

from sys import argv
from parsers import parse_fasta, LongFastaID
from writers import write_fasta

if __name__ == "__main__":
    fasta = argv[1]
    out_file = argv[2]

    fasta_dict = parse_fasta(fasta)
    prot_id_to_seq = {}

    for fasta_id, sequence in fasta_dict.items():
        fasta_id = LongFastaID(fasta_id)
        prot_id = fasta_id.protein_id
        assert prot_id not in prot_id_to_seq
        prot_id_to_seq[prot_id] = sequence

    write_fasta(prot_id_to_seq, out_file)

        
