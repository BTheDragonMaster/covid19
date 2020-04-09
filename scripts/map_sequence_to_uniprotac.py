#!/usr/bin/env python

from parsers import parse_fasta, BlastHit, parse_blast_output
from find_unique_sequences import reverse_fasta_dict
import subprocess
import os

dirname = os.path.dirname(__file__)
CODE_DIR = os.path.join(dirname, '../data/code_to_accesion.fasta')
UNIQUE_SEQ_DIR = os.path.join(dirname, '../data/unique_sequences.fasta')
REFERENCE_PROTEOME = os.path.join(dirname, '../data/reference_proteome_complete.fasta')
TEMP_BLAST = os.path.join(dirname, '../data/temp/temp_blast_results.txt')

def run_muscle(in_file):
    out_file = in_file.split('.')[0] + '_aligned.faa'
    command = ['muscle', '-in', in_file, '-out', out_file]
    subprocess.check_call(command)

def run_blastp(in_file):
    subjects = REFERENCE_PROTEOME
    queries = in_file
    command = ['blastp', '-query', queries, '-subject', subjects, '-out',
               TEMP_BLAST, '-outfmt',
               "6 qseqid sseqid pident length mismatch qstart qend qlen sstart \
send slen evalue bitscore qcovs"]
    subprocess.check_call(command)

    

if __name__ = "__main__":
    fasta = argv[1]

    unique_id_to_seq = parse_fasta(UNIQUE_SEQ_DIR)
    unique_seq_to_id = reverse_fasta_dict(unique_id_to_seq)
    
    new_id_to_seq = parse_fasta(fasta)
    new_seq_to_id = reverse_fasta_dict(new_id_to_seq)

    for sequence in new_seq_to_id:
        if sequence in unique_seq_to_id:
            
            
            

    
    pass



               '
    


    
