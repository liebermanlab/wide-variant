import os, sys, itertools, 
from Bio import SeqIO

def get_sample_names(): 
    # Read in samples.csv file 
    sn = []
    align_types = [] 
    with open('samples.csv') as s: 
        next(s)
        for line in s: 
            entry = line.strip().split(',')
            sn.append(entry[3]) 
            align_types.append(entry[4])
    return sn, align_type 

def get_alignment_filter_map(): 
    align_to_filter = {}
    with open('alignment_params.csv', 'rU') as a: 
        next(a) 
        for line in a: 
            entry = line.strip().split(',') 
            align_name, filter_name = entry[0], entry[1]
            align_to_filter[align_name] = filter_name 
    return align_to_filter

def append_sample_unaligned_reads(sname, salign, afmap, fh1, fh2): 
    fastq1_name = sname+'/'+afmap[salign]+'/unaligned.1.fastq'
    fastq2_name = sname+'/'+afmap[salign]+'/unaligned.2.fastq' 
    
    print 'Adding unaligned reads from %s and %s' %(fastq1_name, fastq2_name) 
    fastq_iter1 = SeqIO.parse(open(fastq1_name), 'fastq') 
    fastq_iter2 = SeqIO.parse(open(fastq2_name), 'fastq') 
    for rec1, rec2 in itertools.izip(fastq_iter1, fastq_iter2):
        SeqIO.write(rec1, fh1, 'fastq')
        SeqIO.write(rec2, fh2, 'fastq') 
        

def main(): 
    # get all sample names
    sample_names, sample_align = get_sample_names()
    # get alignment params and corresponding filter
    align_filter_map = get_alignment_filter_map() 

    # initialize master file
    unaligned_1 = 'unaligned.all.1.fastq'
    unaligned_2 = 'unaligned.all.2.fastq' 
    if os.path.isfile(unaligned_1): 
        os.remove(unaligned_1)
        print '\nRemoved existing file %s' %unaligned_1
    if os.path.isfile(unaligned_2): 
        os.remove(unaligned_2) 
        print '\nRemoved existing file %s' %unaligned_2
    fh1 = open(unaligned_1, 'a') 
    fh2 = open(unaligned_2, 'a') 

    # iterate through all files 
    for i in range(len(samplenames)): 
        append_sample_unaligned_reads(sample_names[i], sample_align[i], align_filter_map, fh1, fh2) 
    print '\nFinished appending all unaligned reads' 
    fh1.close()
    fh2.close()

if __name__ == '__main__': 
    print '\n Merging all unaligned fastq files' 
    main() 
