#!/usr/bin/env python

import subprocess
import os
from sys import argv

def run_muscle(in_file):
    out_file = in_file.split('.')[0] + '_aligned.faa'
    command = ['muscle', '-in', in_file, '-out', out_file]
    subprocess.check_call(command)

if __name__ == "__main__":
    in_dir = argv[1]
    for fasta in os.listdir(in_dir):
        if fasta[-6:] == '.fasta':
            fasta_file = in_dir + fasta
            run_muscle(fasta_file)
    
    
    
