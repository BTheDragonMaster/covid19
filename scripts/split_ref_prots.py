#!/usr/bin/env python

from parsers import parse_fasta, LongFastaID
from writers import write_fasta

from sys import argv


REFACC = "NC_045512.2"

def separate_ref_from_nonref(fasta_dir):
    ref_fasta_dict = {}
    nonref_fasta_dict = {}
    
    fasta_dict = parse_fasta(fasta_dir)
    for seq_id, sequence in fasta_dict.items():
        seq_id = LongFastaID(seq_id)
        if seq_id.genome_acc == REFACC:
            ref_fasta_dict[seq_id.protein_id] = sequence
        else:
            nonref_fasta_dict[seq_id.protein_id] = sequence

    return ref_fasta_dict, nonref_fasta_dict

if __name__ == "__main__":
    fasta_dir = argv[1]

    ref_fasta_dict, nonref_fasta_dict = separate_ref_from_nonref(fasta_dir)
    write_fasta(ref_fasta_dict, 'reference_proteome.fasta')
    write_fasta(nonref_fasta_dict, 'non-reference_covid19_proteins.fasta')
    


            
        
    
