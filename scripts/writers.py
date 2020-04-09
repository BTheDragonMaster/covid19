#!/usr/bin/env python

def write_residue_mapping(mapping, out_dir):
    query_residues = sorted(mapping.keys())
    
    with open(out_dir, 'w') as out_file:
        for res in query_residues:
            ref_res = mapping[res]
            out_file.write(f'{res}\t{ref_res}\n')
        

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


def write_code_to_accession(code_to_accessions, out_dir):
    codes = sorted(list(code_to_accessions.keys()))
    with open(out_dir, 'w') as out_file:
        for code in codes:
            out_file.write(code)
            out_file.write('\t')
            accessions = code_to_accessions[code]
            nr_of_accessions = len(accessions)
            for i, accession in enumerate(accessions):
                if i == len(accessions) - 1:
                    
                    out_file.write(f'{accession}')
                else:
                    out_file.write(f'{accession},')
            out_file.write('\n')
                
        


        
    
