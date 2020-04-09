#!/usr/bin/env python

from operator import itemgetter, attrgetter

from parsers import parse_blast_output, BlastHit, BlastResult

def get_best_hits(query_to_combined_subject):
    query_to_best_hit = {}
    for query, subject_to_combined_hit in query_to_combined_subject.items():
        query_to_best_hit[query] = sort_hits(subject_to_combined_hit)[0]

    return query_to_best_hit
    
    

def sort_hits(subject_to_combined_hit):
    all_hits = []
    for subject, combined_hit in subject_to_combined_hit.items():
        all_hits.append((combined_hit, subject))

    
    sorted_hits = sorted(all_hits, key = lambda x: attrgetter('qcov', 'scov')(x[0]), reverse = True)
#    sorted_hits = sorted(all_hits, key = lambda x: attrgetter('scov', 'qcov')(x[0]), reverse = True)
#    sorted_hits = sorted(sorted_hits, key = lambda x: attrgetter('worst_eval')(x[0]))
    

    return sorted_hits
    
    

def order_hits(hit_dict):
    new_hit_dict = {}
    for query_id, hits in hit_dict.items():
        
        
        subject_to_hits = {}
        for hit in hits:
            subject_id = hit.sseqid
            if not subject_id in subject_to_hits:
                subject_to_hits[subject_id] = []
                
            subject_to_hits[subject_id].append(hit)

        new_hit_dict[query_id] = subject_to_hits
            

    return new_hit_dict

def combine_all_hits(query_to_hit_dict):
    combined_hit_dicts = {}
    
    for query, subject_to_hits in query_to_hit_dict.items():
        combined_hit_dict = combine_hits(subject_to_hits)
        combined_hit_dicts[query] = combined_hit_dict

    return combined_hit_dicts
        


def combine_hits(subject_to_hits):
    subject_to_combined_hit = {}
    for subject_id, hits in subject_to_hits.items():
        subject_to_combined_hit[subject_id] = BlastResult(hits)

    return subject_to_combined_hit

def filter_hits(hit_dict, fasta_dict):
    filtered_hits = {}
    identical_hits = {}
    rejected_hits = {}
    for query_id, subject_to_combined_hit in hit_dict.items():
        filtered_hits[query_id] = {}
        identical_hits[query_id] = {}
        rejected_hits[query_id] = {}
        for subject_id, combined_hit in subject_to_combined_hit.items():
            
            subject_length = len(fasta_dict[subject_id])
            subject_aa_coverage = float(combined_hit.cover) / subject_length
            subject_length_coverage = float(combined_hit.length) / subject_length

            smallest_ident = combined_hit.get_smallest_ident()
            smallest_eval = combined_hit.get_smallest_eval()

            if smallest_ident < 100 and subject_aa_coverage > 0.8 and subject_length_coverage > 0.8 and smallest_eval < 0.001 :
                filtered_hits[query_id][subject_id] = combined_hit
            elif int(smallest_ident) == 100 and smallest_eval < 0.001:
                identical_hits[query_id][subject_id] = combined_hit
            else:
                rejected_hits[query_id][subject_id] = combined_hit

    return filtered_hits, identical_hits, rejected_hits
