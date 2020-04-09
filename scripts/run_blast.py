import os
from sys import argv
import subprocess

dirname = os.path.dirname(__file__)

REFERENCE_PROTEOME = os.path.join(dirname, '../data/reference_proteome_complete.fasta')

def run_blastp(in_file, out_file):
    subjects = REFERENCE_PROTEOME
    queries = in_file
    command = ['blastp', '-query', queries, '-subject', subjects, '-out',
               out_file, '-outfmt',
               "6 qseqid sseqid pident length mismatch qstart qend qlen sstart \
send slen evalue bitscore qcovs"]
    subprocess.check_call(command)

if __name__ == "__main__":
    run_blastp(argv[1], argv[2])
