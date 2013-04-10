'''
HC 4/9/2013
Takes a VCF file and edits a reference genome (fasta or genbank) to replace SNPs. 
Write output in new fasta or genbank file. 
'''

import sys, os 
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def load_genome_sequence(filename, filetype):
    # load gb file 
    seq_record = SeqIO.read( filename, filetype )
    print "Loaded genome %s, which is %i bp long" %(seq_record.id, len(seq_record))
    # convert Seq file to mutable
    mutable_seq = seq_record.seq.tomutable() 
    return seq_record, mutable_seq

def edit_genome_sequence( mutable_seq, var_file ):
    # params
    nucleotides = ('A', 'T', 'C', 'G')
    counter = 0

    # open variant file
    var_filehandle = open(variant_file, "rU")

    # for each variant, replace reference genome seq 
    for line in var_filehandle: 
        pos, ref, var, qual = line.split()
        if float(qual)>30: # qual filter 
            if (ref and var) in nucleotides: # replace only SNPs, not indels
                fasta_nucleotide = mutable_seq[int(pos)-1] 
                # check that ref in VCF and FASTA are the same 
                if ref == fasta_nucleotide: 
                    mutable_seq[int(pos)-1] = var 
                    counter +=1 
    var_filehandle.close()

    print 'Replaced %i SNP positions' %counter
    return mutable_seq

def save_new_file( edited_mutseq, unedited_seq_record, variant_file, filetype, output_file ):
    # save edited Seq as gb or fasta file 
    edited_seq = edited_mutseq.toseq() 
    # make new seq record
    new_description = unedited_seq_record.description+' Edited with '+variant_file
    edited_seq_record = SeqRecord( 
            edited_seq,                                     # mutated sequence
            id=unedited_seq_record.id,                      # id 
            features=unedited_seq_record.features,          # gene features 
            annotations=unedited_seq_record.annotations,    # annotations
            description=new_description                     # description; modified 
            )

    # write new genome from new seq record 
    out_handle = open(output_file, 'w')
    SeqIO.write(edited_seq_record, out_handle, filetype) 
    out_handle.close() 

if __name__ == '__main__': 
    # ___ obtain command line arguments ___ #
    if len(sys.argv) != 4: 
        print len(sys.argv)
        print 'Usage: \n\t python edit_reference_with_vcf.py [fasta file | gb file] [vcf file] [new fasta or gb file name]'
        sys.exit(1) 
    else:  
        # get genome file 
        genome_file = sys.argv[1]
        if genome_file.endswith('.fasta'):
            filetype = "fasta"
        elif genome_file.endswith(('.gb', 'gbk')):
            filetype = "genbank"
        else:
            print 'Invalid genome file type. Please try again.' 
            sys.exit(1)
        # vcf file 
        variant_file = sys.argv[2] 
        # new output file name 
        output_file = sys.argv[3]
    
    # get genome sequence as mutable string Seq object
    old_seq_record, mut_seq = load_genome_sequence( genome_file, filetype ) 
    
    # edit each position 
    edited_mutseq = edit_genome_sequence( mut_seq, variant_file ) 
    
    # make new file 
    save_new_file( edited_mutseq, old_seq_record, variant_file, filetype, output_file ) 
