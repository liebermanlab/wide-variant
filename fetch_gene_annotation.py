import os, sys
from Bio import Entrez, SeqIO

def fetch_gene(geneid):
    record = [] 
    h = Entrez.efetch(db="gene", id=str(geneid), rettype="gb", retmode="text") 
    for line in h: record.append(line.strip())
    annotation = record[2] # HARD CODING 
    #print '\nGene %s, annotation %s\t' %(geneid, annotation) 
    if annotation.find('['):
        annotation_return = annotation[0:annotation.find('[')]
    else: 
        annotation_return = annotation 
    return annotation_return 

if __name__ == '__main__': 
    Entrez.email = 'hattiechung@gmail.com' 
    gid = sys.argv[1] 
    ann = fetch_gene(gid) 
    sys.stdout.write(ann) 
