#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 28 18:51:07 2020

@author: tamilieberman
"""

'''
straincraker.py reassigns reads from a kraken output file and produces a new kraken output file

It uses a purity threshold to reassign taxa to the taxonomic level at which the percent
of kmers assigned at a node or its descendents is above the specified the specified treshold
Kmers above the taxonomic level in consideration are not included in this calculation
Requires a precompiled tree, generated from straincracker_loadtree.py and nodes.dmp

Copyright (C) 2020 Tami Lieberman, tami@mit.edu


Required Parameters:
   -i,--infile X....................kraken output file, gzipped
   -t,--treefile X.....................straincraker treefile
   -o, --output X......................modified kraken report file
   -p, --purity X......................purity threshold 
   
Example: python straincracker.py --t treematrix.npy --i Meta-F_Extract-Single_Ch_krakSeq.txt.gz --o Meta-F_Extract-Single_Ch_krakSeq.txt.cracked --p .99


#steps to all of straincraker
#make treematrix, krakendatabase format
#kraken
#strainkrak
#convert
#bracken
'''

# %% Run the actual program




#################################################################################

import numpy as np
import gzip
import argparse
import os
#################################################################################

# %%
def main():
    #Parse arguments
    parser = argparse.ArgumentParser(description='Reassign reads using treefile from strainkraker_loadtree.py')
    parser.add_argument('--i', metavar='infile', nargs='?', required=True, help='gzipped output from kraken')
    parser.add_argument('--o', metavar='outfile', nargs='?', required=True, help='destination output file')
    parser.add_argument('--t', metavar='treefile', nargs='?', required=True,  help='source tree file')
    parser.add_argument('--p', metavar='purity', type=float, required=True, nargs='?', help='purity threshold')
    args = parser.parse_args()
    strainCrack(args.t,args.i,args.o,args.p)
    
    
#%%

def read_in_results(kstrings):
    
    kcounts=dict({0:0, 1:0})
        
    for i,s in enumerate(kstrings):
        x=s.split(':')
        if x[0] != '|':
            if int(x[0]) in kcounts:
                kcounts[int(x[0])]+=int(x[1])
            else:
                kcounts[int(x[0])]=int(x[1])
    
    kcounts.pop(0) #ignore kmers alinging to levels 0 or 1
    kcounts.pop(1) #ignore kmers alinging to levels 0 or 1

    t=np.array(list(kcounts.keys()))
    k=np.array(list(kcounts.values()))
    
    return t,k


# %%

def classify_read(taxa,kmers, treem, purity_thresh,tlevels):
    
    found_taxa=taxa.copy()
    kmer_values=kmers.copy()
    taxa_categorized=0
    subtree=treem[found_taxa,:];
    
    num_levels=max(np.where(subtree>0)[1])

    #descend tree
    #calculate kmers at each level  
    for l in range(1,num_levels+1):

        taxa_at_this_level_or_below=np.logical_and(kmer_values>0,tlevels>=l)
        
        #first try to classify at this level
        if np.sum(taxa_at_this_level_or_below) == 1:
            taxa_categorized=found_taxa[taxa_at_this_level_or_below][0]
            break
        elif np.sum(taxa_at_this_level_or_below) < 1:
            #go with previous assignment
            break
        else: 
         #need to coninue classifying

            #find all taxa at this level that have descendents w kmers
            level_i_classification_for_taxa_found=subtree[:,l]
            taxa_at_this_level_w_evidence=np.unique(level_i_classification_for_taxa_found[level_i_classification_for_taxa_found>0])    
                    
            #are there 0, 1, or more taxa at this level?
            if taxa_at_this_level_w_evidence.size == 1:
                #just 1 taxa, pure by definitaion, continue going down tree
                taxa_categorized = taxa_at_this_level_w_evidence[0]
            elif taxa_at_this_level_w_evidence.size < 1:
                #report prev taxa
                print('Warning: somehow got here')
                break
            else: # taxa_at_this_level_w_evidence.size > 1:
                
                #check which is best path to go down, or report prev
                kmers_at_level=np.zeros(np.shape(taxa_at_this_level_w_evidence),dtype=int)
                for i,t in enumerate(taxa_at_this_level_w_evidence):
                    self_and_children_taxa=np.equal(level_i_classification_for_taxa_found,t);
                    kmers_at_level[i]=sum(kmer_values[self_and_children_taxa])
    
                purity=max(kmers_at_level)/sum(kmers_at_level)
                if purity > purity_thresh:                
                    
                    taxa_categorized=taxa_at_this_level_w_evidence[kmers_at_level==max(kmers_at_level)]
                    taxa_categorized = taxa_categorized[0]
                    
                    #set all counts on paths not below this to 0
                    to_delete=np.not_equal(level_i_classification_for_taxa_found,taxa_categorized)
                    kmer_values[to_delete]=0

                else:
                    #report prev taxa
                    break
        
    
    return taxa_categorized



# %%
def strainCrack(treefile,infile,outfile,purity_threshold):
    
    #i mport tree structure
    tree = np.load(treefile)
    
        
    # make an array that says what taxonmic level each taxa is at, helpful for later
    taxa_nums=np.array([range(tree.shape[0])])
    taxa_exists=np.sum(np.equal(tree,taxa_nums.T)>0,1)
    taxonomic_levels=np.zeros(taxa_exists.shape,dtype=int)
    taxonomic_levels[taxa_exists>0]=np.array((np.nonzero(np.equal(tree,taxa_nums.T))))[1]
    
    
    # read in input and output file

    f=gzip.open(infile,"rt")
    of=open(outfile,"w")
    l=f.readline()
    
    # test
    changed=0
    tested=0
    while len(l) > 0 and l != "\n":
        line=l.strip().split()
        
        if line[0]=='C':
            taxa, kmers =read_in_results(line[4:])
            if kmers.size > 1: 
                tested+=1
                new_classification=classify_read(taxa,kmers,tree,purity_threshold, taxonomic_levels[taxa])
                if int(line[2]) != new_classification:
                    changed+=1
                    line[2]=str(new_classification)
            newline1='\t'.join(map(str, line[:4]))
            newline2=' '.join(map(str, line[4:]))
            of.write(newline1 + '\t' + newline2 + '\n')         
        else:
            of.write(l)
            
        l=f.readline()
 
    of.close()
    f.close()
        
    return changed, tested

if __name__ == "__main__":
    main()