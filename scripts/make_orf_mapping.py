#!/usr/bin/env python

from parsers import parse_fasta, LongFastaID
from sys import argv

if __name__ == "__main__":
    fasta = argv[1]
    orf_dir = argv[2]
    
    fasta_dict = parse_fasta(fasta)
    with open(orf_dir, 'w') as orf_file:
        for seq_id in fasta_dict:
            seq_id = LongFastaID(seq_id)
            prot_id = seq_id.protein_id
            orf = seq_id.protein
            orf_file.write(f'{prot_id}\t{orf}\n')
        
    
