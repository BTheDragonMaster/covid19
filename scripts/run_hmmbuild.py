#!/usr/bin/env python

import subprocess
import os
from sys import argv

def run_buildhmm(in_file):
    out_file = in_file.split('.')[0] + '.hmm'
    command = ['hmmbuild', out_file, in_file]
    subprocess.check_call(command)

if __name__ == "__main__":
    in_dir = argv[1]
    for msa in os.listdir(in_dir):
        if msa[-4:] == '.faa':
            msa_file = in_dir + msa
            run_buildhmm(msa_file)
