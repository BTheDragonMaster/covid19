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

class BlastResult:
    def __init__(self, hits):
        self.hits = hits
        self.set_subject_ranges()
        self.set_scov()
        self.set_qcov()
        self.set_length()
        self.set_worst_eval()
        self.set_smallest_ident()
        print(self.subject_ranges)

    def set_subject_ranges(self):
        subject_starts = []
        subject_ends = []
        self.subject_ranges = []
        for hit in self.hits:
            subject_starts.append(hit.sstart)
            subject_ends.append(hit.send)
            self.subject_ranges.append((hit.sstart, hit.send))

        self.subject_start = min(subject_starts)
        self.subject_end = max(subject_ends)
        self.subject_ranges.sort(key = lambda x: x[0])
        self.filter_overlapping_ranges()

    def filter_overlapping_ranges(self):
        overlap = self.subject_ranges[0]
        overlaps = []
        counter = 0
        for start, end in self.subject_ranges:
            if counter > 1:
                if overlap[1] >= start >= overlap[0]:
                    overlap = (overlap[0], max(end, overlap[1]))
                else:
                    overlaps.append(overlap)
                    overlap = (start, end)
            else:
                pass
            counter += 1
        overlaps.append(overlap)

        self.subject_ranges = overlaps
            
    def set_scov(self):
        scov = 0
        for start, end in self.subject_ranges:
            cover = end - start + 1
            scov += cover

        self.scov = 100 * scov / float(self.hits[0].slen)

    def set_qcov(self):
        self.qcov = self.hits[0].qcov
    
    def set_length(self):
        self.length = self.subject_end - self.subject_start + 1

    def set_smallest_ident(self):
        identities = []
        for hit in self.hits:
            identities.append(hit.pident)

        self.smallest_ident = min(identities)

    def set_worst_eval(self):
        evals = []
        for hit in self.hits:
            evals.append(hit.evalue)

        self.worst_eval  = max(evals)

class BlastHit:

    def __init__(self, line):
        self.parse_line(line)

    def parse_line(self, line):
        self.qseqid, self.sseqid, self.pident, self.length,\
                     self.mismatch, self.qstart,\
                     self.qend, self.qlen, self.sstart, self.send, \
                     self.slen, self.evalue,\
                     self.bitscore, self.qcov = line.split('\t')

        self.pident = float(self.pident)
        self.length = int(self.length)
        self.mismatch = int(self.mismatch)
        self.qlen = int(self.qlen)
        self.qstart = int(self.qstart)
        self.qend = int(self.qend)
        self.sstart = int(self.sstart)
        self.send = int(self.send)
        self.slen = int(self.slen)
        self.evalue = float(self.evalue)
        self.bitscore = float(self.bitscore)
        self.qcov = float(self.qcov)
        self.scov = self.length / float(self.slen)

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
        


        
    
