#!/usr/bin/env python

def parse_mapping(mapping_dir):
    placeholder_to_id = {}
    with open(mapping_dir, 'r') as mapping_file:
        for line in mapping_file:
            line = line.strip()
            placeholder, id = line.split('\t')
            placeholder_to_id[placeholder] = id

    return placeholder_to_id

def parse_blast_output(blast_dir):
    hit_dict = {}
    with open(blast_dir, 'r') as blast_file:
        for line in blast_file:
            line = line.strip()
            blast_hit = BlastHit(line)
            if not blast_hit.qseqid in hit_dict:
                
                hit_dict[blast_hit.qseqid] = []

            hit_dict[blast_hit.qseqid].append(blast_hit)

    return hit_dict

class BlastHit:

    def __init__(self, line):
        self.parse_line(line)

    def parse_line(self, line):
        self.qseqid, self.sseqid, self.pident, self.length,\
                     self.mismatch, self.gapopen, self.qstart,\
                     self.qend, self.sstart, self.send, self.evalue,\
                     self.bitscore = line.split('\t')

        self.pident = float(self.pident)
        self.length = int(self.length)
        self.mismatch = int(self.mismatch)
        self.gapopen = int(self.gapopen)
        self.qstart = int(self.qstart)
        self.qend = int(self.qend)
        self.sstart = int(self.sstart)
        self.send = int(self.send)
        self.evalue = float(self.evalue)
        self.bitscore = float(self.bitscore)
        self.scover = self.length

def parse_id_from_prot_file(prot_id):
    prot_id = prot_id.split('|')[0].strip()
    return prot_id
    
def parse_fasta(fasta_dir):
    with open(fasta_dir, 'r') as fasta_file:
        fasta_dict = {}
        sequence = []
        ID = None
        for line in fasta_file:
            line = line.strip()
            if line.startswith('>'):
                if ID:
                    sequence = ''.join(sequence)
                    fasta_dict[ID] = sequence
                    sequence = []
                    
                ID = line[1:]
                
            else:
                sequence.append(line)
        sequence = ''.join(sequence)
        fasta_dict[ID] = sequence

    return fasta_dict

class LongFastaID:

    def __init__(self, fasta_id):
        self.gene = None
        self.locus_tag = None
        self.db_xref = None
        self.protein = None
        self.exception = None
        self.protein_id = None
        self.location = None
        self.gbkey = None
        self.partial = None
        self.frame = None
        self.transl_except = None
        
        self.parse_fasta_ID(fasta_id)

    def parse_fasta_ID(self, ID):
        accessions = ID.split()[0].split('|')[1]
        self.genome_acc, self.prot_acc = accessions.split('_prot_')

        other_info = ID.split(' [')[1:]
        for i, data_point in enumerate(other_info):
            data_type, data_value = data_point[:-1].split('=')
            if data_type == 'gene':
                self.gene = data_value
            elif data_type == 'locus_tag':
                self.locus_tag = data_value
            elif data_type =='protein':
                self.protein = data_value.replace(' ', '-')
            elif data_type == 'db_xref':
                self.db_xref = data_value
            elif data_type == 'exception':
                self.exception = data_value
            elif data_type == 'protein_id':
                self.protein_id = data_value
            elif data_type == 'location':
                self.location = data_value
            elif data_type == 'gbkey':
                self.gbkey = data_value
            elif data_type == 'frame':
                self.frame = data_value
            elif data_type == 'partial':
                self.partial = data_value
            elif data_type == 'transl_except':
                self.transl_except = data_value
            else:
                print(data_type)
        


        
    
