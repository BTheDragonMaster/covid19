#!/usr/bin/env python

from sys import argv
from parsers import parse_blast_output, BlastHit, parse_fasta_simple, FastaID, ProteinRecord, parse_mapping
from writers import write_fasta

def sort_hits_by_subject(hit_dicts):
    subject_to_hits = {}
    for hit_dict in hit_dicts:
        for id, hits in hit_dict.items():
            for hit in hits:
                subject = hit.sseqid
                if not subject in subject_to_hits:
                    subject_to_hits[subject] = []
                subject_to_hits.append(hit)

    return subject_to_hits

def find_polyprotein_fragments(subject_to_hits):
    for subject, hits in subject_to_hits.items():
        rangers = []
        for hit in hits:
            query_id = 
        


if __name__ == "__main__":
    blast_output_nonref = argv[1]
    blast_output_pdb = argv[2]
    
    reference_proteome = argv[3]
    
    pdb_fasta = argv[5]
    nonref_proteins = argv[6]
    
    nonref_mapping = argv[7]

    nonref_mapping = parse_mapping(nonref_proteins)
    ref_mapping = parse_mapping(reference_mapping)

    hit_dicts = []

    hit_dict_pdb = parse_blast_output(blast_output_pdb, subject_mapping=ref_mapping)
    hit_dict_nonref = parse_blast_output(blast_output_nonref, query_mapping=nonref_mapping,
                                         subject_mapping=ref_mapping)

    hit_dicts.append(hit_dict_pdb)
    hit_dicts.append(hit_dict_nonref)

    subject_to_hits = sort_hits_by_subject(hit_dicts)

    
    
