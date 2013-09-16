from Bio import SeqIO
import os, sys, itertools

def main(ref, contig, output): 
    # initialize output 
    output_fh = open(output, 'w')
    # add reference entries
    ref_fh = SeqIO.parse(open(ref), 'fasta')
    for rec1 in ref_fh:
        SeqIO.write(rec1, output_fh, 'fasta')
    # add contig entries
    contig_fh = SeqIO.parse(open(contig), 'fasta') 
    contig_len = 0
    for rec2 in contig_fh: 
        SeqIO.write(rec2, output_fh, 'fasta') 
        contig_len += len(rec2) 
    # close handle 
    output_fh.close() 

    print '\nContigs have total %i characters' %contig_len 

if __name__ == '__main__': 
    args = sys.argv[1:] 
    if len(args) is not 3: 
        print '\n\nUsage: python add_contig_to_reference.py [Reference fasta file] [contig fasta file] [output name]'
        sys.exit(1) 
    [ref_file, contig_file, output_file] = args

    main(ref_file, contig_file, output_file) 
