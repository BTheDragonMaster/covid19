import xml.etree.cElementTree as ET
import requests as rq

def get_corona_nuc_ids():
    aElement = ET.fromstring(rq.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nucleotide&term=%22Severe%20acute%20respiratory%20syndrome%20coronavirus%202%22&RetMax=1000000").text)
    return [x.text for x in aElement.iter("Id")]

def get_corona_prot_ids():
    aElement = ET.fromstring(rq.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=protein&term=%22Severe%20acute%20respiratory%20syndrome%20coronavirus%202%22&RetMax=1000000").text)
    return [x.text for x in aElement.iter("Id")]
    
def write_covid_nuc_sequences(fname):
    ids = get_corona_nuc_ids()
    fasta = rq.get(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id={','.join(ids)}&rettype=fasta").text
    f = open(fname, "w")
    f.write(fasta)
    f.close()

def write_covid_prot_sequences(fname):
    ids = get_corona_prot_ids()
    fasta = list()
    for i in range(0, len(ids), 20):
        fasta.append(rq.get(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id={','.join(ids[i:i+20])}&rettype=fasta").text)
    f = open(fname, "w")
    f.write("/n".join(fasta))
    f.close()

if __name__ == "__main__":
    write_covid_nuc_sequences("covid_seqs_nuc.fasta")
    write_covid_prot_sequences("covid_seqs_prot.fasta")


