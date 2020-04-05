from parsers import parse_mapping, parse_blast_output, parse_fasta_simple, FastaID
from sys import argv

def make_query_mapping(query_fasta):
    mapping_dict = {}
    fasta_dict = parse_fasta_simple(query_fasta)
    for ID in fasta_dict:
        new_ID = ID.split(':')[0]
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

def make_protid_mapping(subject_mapping):
    protid_mapping = {}
    for placeholder, ID in subject_mapping.items():
        fastaID = FastaID(ID)
        protid_mapping[ID] = fastaID.protein_id

    return protid_mapping
    

if __name__ == "__main__":
    blast_output = argv[1]
    ref_mapping = argv[2]
    pdb_fasta = argv[3]

    query_mapping = make_query_mapping(pdb_fasta)
    subject_mapping = parse_mapping(ref_mapping)

    
    
    hit_dict = parse_blast_output(blast_output, query_mapping, subject_mapping)
    protid_mapping = make_protid_mapping(subject_mapping)
    print(protid_mapping)
    write_hits(hit_dict, protid_mapping, 'pdb_to_refseqs.txt')


