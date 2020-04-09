from parsers import parse_mapping, parse_blast_output, parse_fasta, FastaID
from sys import argv

import os

dirname = os.path.dirname(__file__)

REFERENCE_PROTEOME = os.path.join(dirname, '../data/reference_proteome_complete.fasta')
PDB_SEQUENCES = os.path.join(dirname, '../data/pdb_sequences_full.fasta')

def make_query_mapping(query_fasta):
    mapping_dict = {}
    fasta_dict = parse_fasta(query_fasta)
    for ID in fasta_dict:
        new_ID = ID.split('|')[0]
        print(new_ID)
        mapping_dict[ID] = new_ID

    return mapping_dict

def write_hits(hit_dict, subject_mapping, out_dir):
    with open(out_dir, 'w') as out_file:
        for ID, hits in hit_dict.items():
            out_file.write(ID)
            sseqids = set([])
            for hit in hits:
                if hit.evalue < 0.001:
                    sseqids.add(subject_mapping[hit.sseqid])

            for sseqid in sseqids:

                out_file.write(f'\t{sseqid}')
            out_file.write('\n')
    

if __name__ == "__main__":
    blast_output = argv[1]
    
    hit_dict = parse_blast_output(blast_output)
    
    print(protid_mapping)
    write_hits(hit_dict, protid_mapping, 'pdb_to_refseqs.txt')


