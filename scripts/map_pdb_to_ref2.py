#!/usr/bin/env python

from parsers import parse_fasta, parse_blast_output, BlastHit, BlastResult
from handle_blastp_results import get_best_hits, order_hits, sort_hits, combine_all_hits, combine_hits, filter_hits
from run_muscle import run_muscle
from writers import write_fasta, write_residue_mapping

import os


dirname = os.path.dirname(__file__)

BLAST_OUTPUT = os.path.join(dirname, '../data/pdb_vs_ref_blastp.txt')
MAPPING_DIR = os.path.join(dirname, '../data/mappings/pdb_to_ref/')
PDB_SEQS_FULL = os.path.join(dirname, '../data/pdb_sequences_full.fasta')
PDB_SEQS_STRUCTURE = os.path.join(dirname, '../data/pdb_sequences_structure.fasta')
REFERENCE_PROTEOME = os.path.join(dirname, '../data/reference_proteome_complete.fasta')
TEMP = os.path.join(dirname, '../data/temp/muscle/')

def get_query_and_ref_ids(aligned_fasta_dict):
    query_id = None
    ref_id = None
    for fasta_id in aligned_fasta_dict:
        if not fasta_id.startswith('YP') and not fasta_id.startswith('P0'):
            query_id = fasta_id
        else:
            ref_id = fasta_id

    

    return query_id, ref_id
    

def make_mapping(aligned_fasta_dict):
    query_id, ref_id = get_query_and_ref_ids(aligned_fasta_dict)

    query_seq = aligned_fasta_dict[query_id]
    ref_seq = aligned_fasta_dict[ref_id]

    mapping = {}
    
    query_res = 0
    ref_res = 0
    
    for i, aa in enumerate(query_seq):
        ref_aa = ref_seq[i]
        if aa != '-':
            query_res += 1
        

        if ref_aa != '-':
            ref_res += 1

        if aa != '-' and ref_aa != '-':
            mapping[query_res] = ref_res

    return mapping
            

if __name__ == "__main__":

    blast_results = parse_blast_output(BLAST_OUTPUT)

    id_to_sequence = parse_fasta(PDB_SEQS_STRUCTURE)
    ref_to_sequence = parse_fasta(REFERENCE_PROTEOME)

    ordered_blast_results = order_hits(blast_results)
    combined_blast_results = combine_all_hits(ordered_blast_results)
    best_hits = get_best_hits(combined_blast_results)

    for query, best_hit in best_hits.items():
        fasta_dict = {}
        ref_id = best_hit[1]
        fasta_dict[query] = id_to_sequence[query]
        fasta_dict[ref_id] = ref_to_sequence[ref_id]
        
        temp_fasta = f'{TEMP}{query}.fasta'
        write_fasta(fasta_dict, temp_fasta)
        temp_aligned = f'{TEMP}{query}_aligned.fasta'
 #       if not os.path.isfile(temp_aligned):
        run_muscle(temp_fasta, temp_aligned)
        aligned_fasta_dict = parse_fasta(temp_aligned)
        mapping = make_mapping(aligned_fasta_dict)
        mapping_out = f'{MAPPING_DIR}{query}_{ref_id}.txt'
        write_residue_mapping(mapping, mapping_out)
        
        
        
        
 #       print(query, best_hit[1], best_hit[0].qcov, best_hit[0].scov, best_hit[0].worst_eval, best_hit[0].smallest_ident)
    

    
    

    
    
