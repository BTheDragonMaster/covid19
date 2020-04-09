#!/usr/bin/env python

from sys import argv
from parsers import parse_blast_output, BlastHit, BlastResult, parse_fasta parse_mapping
from handle_blastp_results import order_hits, combine_hits, filter_hits
from writers import write_fasta
import copy


def map_subject_to_queries(filtered_hits):
    subject_to_queries = {}

    for query_id, subject_to_combined_hit in filtered_hits.items():
        for subject_id, combined_hit in subject_to_combined_hit.items():
            if not subject_id in subject_to_queries:
                subject_to_queries[subject_id] = []

            subject_to_queries[subject_id].append(query_id)

    return subject_to_queries


def make_fasta_dicts(subject_to_queries, covid19_fasta_dict, orf_mapping):
    fasta_dicts = {}

    for subject_id, query_ids in subject_to_queries.items():
        subject_sequence = covid19_fasta_dict[subject_id]
        
        file_name = orf_mapping[subject_id]
        
        new_fasta_dict = {}
        new_fasta_dict[subject_id] = subject_sequence
        
        for query_id in query_ids:
            query_sequence = covid19_fasta_dict[query_id]
            new_fasta_dict[query_id] = query_sequence

        fasta_dicts[file_name] = new_fasta_dict

    return fasta_dicts

def make_seq_to_id_dict(fasta_dict, identical_hits, covid19_fasta_dict):
    seq_to_id = {}
    for seq_id, sequence in fasta_dict.items():
        if not sequence in seq_to_id:
            seq_to_id[sequence] = set([])
        seq_to_id[sequence].add(seq_id)

    for query_id, subject_to_combined_hit in identical_hits.items():
        for subject_id, combined_hit in subject_to_combined_hit.items():
            if subject_id in fasta_dict:
                sequence = covid19_fasta_dict[subject_id]
                seq_to_id[sequence].add(query_id)

    for sequence, ids in seq_to_id.items():
        seq_to_id[sequence] = '_'.join(list(ids))

    return seq_to_id

def make_new_fasta_dict(seq_to_id):
    new_fasta_dict = {}
    for seq, seq_id in seq_to_id.items():
        new_fasta_dict[seq_id] = seq

    return new_fasta_dict
    
            
            
if __name__ == "__main__":
    blast_output = argv[1]
    covid19_fasta = argv[2]
    orf_mapping = argv[3]

    covid19_fasta_dict = parse_fasta(covid19_fasta)
    orf_mapping = parse_mapping(orf_mapping)

    queryid_to_hits = parse_blast_output(blast_output)
    sorted_hit_dict = order_hits(queryid_to_hits)
    combined_hit_dicts = {}
    for query, subject_to_hits in sorted_hit_dict.items():
        combined_hit_dict = combine_hits(subject_to_hits)
        combined_hit_dicts[query] = combined_hit_dict

    filtered_hits, identical_hits, rejected_hits = filter_hits(combined_hit_dicts, covid19_fasta_dict)
 #   print(identical_hits.items())
    
    subject_to_queries = map_subject_to_queries(filtered_hits)

    fasta_dicts = make_fasta_dicts(subject_to_queries, covid19_fasta_dict, orf_mapping)

    for protein_name, fasta_dict in fasta_dicts.items():
        file_name = protein_name + '.fasta'
        seq_to_id = make_seq_to_id_dict(fasta_dict, identical_hits, covid19_fasta_dict)
        id_to_seq = make_new_fasta_dict(seq_to_id)
        write_fasta(id_to_seq, file_name)

    
