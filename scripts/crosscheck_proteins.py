#!/usr/bin/env python

from parsers import parse_fasta, parse_fasta_simple, parse_id_from_prot_file
from sys import argv

def extract_protids_simple(fasta_dict):
    prot_ids = set([])

    for id in fasta_dict:
        prot_ids.add(parse_id_from_prot_file(id))

    return prot_ids

def extract_protids(fasta_dict):
    prot_ids = set([])

    for id, protdata in fasta_dict.items():
        prot_ids.add(protdata.protein_id.split('.')[0])

    return prot_ids
    

if __name__ == "__main__":
    used_fasta = argv[1]
    other_fasta = argv[2]

    fasta_dict_1 = parse_fasta(used_fasta)
    fasta_dict_2 = parse_fasta_simple(other_fasta)

    prot_ids_1 = extract_protids(fasta_dict_1)
    prot_ids_2 = extract_protids_simple(fasta_dict_2)

    print("Extra in used", prot_ids_1 - prot_ids_2)
    for ID in (prot_ids_2 - prot_ids_1):
        print(ID)
