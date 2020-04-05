#!/usr/bin/env python

def write_fasta_from_protref_dict(protref_dict, out_dir):
    seq_id_dir = out_dir.split('.')[0] + '_id_mapping.txt'
    seq_id_file = open(seq_id_dir, 'w')
    with open(out_dir, 'w') as out_file:
        for i, (ID, protref) in enumerate(protref_dict.items()):
            seq_id = i + 1
            out_file.write(f'>{seq_id}\n{protref.sequence}\n')
            seq_id_file.write(f'{seq_id}\t{ID}\n')

    seq_id_file.close()


def write_fasta_from_protref_dict_clean(protref_dict, out_dir):
    with open(out_dir, 'w') as out_file:
        for i, (ID, protref) in enumerate(protref_dict.items()):
            out_file.write(f'>{ID}\n{protref.sequence}\n')

def write_fasta(fasta_dict, out_dir):
    with open(out_dir, 'w') as out_file:
        for id, seq in fasta_dict.items():
            out_file.write(f'>{id}\n{seq}\n')


        
    
