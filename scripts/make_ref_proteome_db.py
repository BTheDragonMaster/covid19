#!/usr/bin/env python
from parsers import parse_mapping, parse_fasta, LongFastaID
from writers import write_fasta
from sys import argv

def get_refseqs(refseq_to_uniprot):
    refseqs = set([])
    for refseq in refseq_to_uniprot:
        refseqs.add(refseq.split('.')[0])

    return refseqs

if __name__ == "__main__":
    fasta = argv[1]
    refseqs = argv[2]

    refseq_to_uniprot = parse_mapping(refseqs)
    refseqs = get_refseqs(refseq_to_uniprot)

    fasta_dict = parse_fasta(fasta)
    refseq_to_seq = {}

    for fasta_id, sequence in fasta_dict.items():
        fasta_id = fasta_id.split('|')[0]
        print(fasta_id)
        fasta_id = fasta_id.strip()
        if fasta_id in refseqs:
            refseq_to_seq[fasta_id] = sequence

    write_fasta(refseq_to_seq, 'reference_proteome_complete.fasta')

    
            
            

    
