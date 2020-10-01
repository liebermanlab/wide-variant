# -*- coding: utf-8 -*-
"""
---Gathers everything together for candidate_mutation_table---

Output:
# path_candidate_mutation_table: where to write
# candidate_mutation_table.mat, ex. results/candidate_mutation_table.mat

---


# Inputs:
     path_to_p_file: where to find all_positions.mat
     path_to_sample_names_file: where to find text file with sample names
         (space delimited)
     path_to_outgroup_boolean_file: where to find text file with outgroup
         booleans (space delimited, 1=outgroup, 0=not)
    path_to_list_of_quals_files: where to find text file with list of
       quals.mat files for each sample (space delimited)
     path_to_list_of_diversity_files: where to find text file with list of
       diversity.mat files for each sample (space delimited)
# Output:
     path_candidate_mutation_table: where to write
     candidate_mutation_table.mat, ex. results/candidate_mutation_table.mat

# Note: All paths should be relative to pwd!


## Version history

     This is adapted from TDL's build_mutation_table_master_smaller_file_size_backup.m
  #   Arolyn, 2018.12.19: This script was written as part of the transition to snakemake. 
          It performs the part of the case step that gathers data for
         Quals and counts and saves candidate_mutation_table.mat
  #   Arolyn, 2019.02.12: Added another matlab variable that stores indel
          statistics called 'indel_counter'.
  #   Tami, 2019.12.12: Converted into python and also added ability save coverage data
"""

print('Importing stuff...')

import numpy as np
import pickle
import scipy.io as sio
import os
import sys
import gzip
from scipy import sparse

def main(path_to_p_file, path_to_sample_names_file, path_to_outgroup_boolean_file, path_to_list_of_quals_files, path_to_list_of_diversity_files, path_to_candidate_mutation_table, path_cov_raw_sparse_matrix,path_cov_norm_sparse_scale_matrix):
    
    print('Starting main...')
    
    pwd=os.getcwd()
    
    # p: positions on genome that are candidate SNPs
    print('Processing candidate SNP positions...')
    
    
    infile=sio.loadmat(path_to_p_file) # from previous step, should include variable called p
    p=infile['p'].flatten()
    p=p-1 #since converting from MATLAB!!!
    print('Total number of positions: ' + str(len(p)))
    
    
    # SampleNames: list of names of all samples
    print('Processing sample names...')
    
    fname =  pwd + '/' + path_to_sample_names_file  
    fid = open( fname, "r" ) # Input is space separated text file, in one line
    SampleNames = fid.readline().split()
    fid.close()
    
    numSamples = len(SampleNames) # save number of samples
    print('Total number of samples: ' + str(numSamples))
    
    
    ## in_outgroup: booleans for whether or not each sample is in the outgroup
    print('Processing outgroup booleans...')
    
    fname =  pwd + '/' + path_to_outgroup_boolean_file  
    fid = open( fname, "r" ) # Input is space separated text file, in one line
    in_outgroup_string = fid.readline().split()
    in_outgroup=np.array(in_outgroup_string)
    fid.close()
    
    
    ## Quals: quality score (relating to sample purity) at each position for all samples
    print('Gathering quality scores at each candidate position...')
    # Import list of directories for where to quals.mat for each sample
    fname =  pwd + '/' + path_to_list_of_quals_files  
    fid = open( fname, "r" ) # Input is space separated text file, in one line
    paths_to_quals_files = fid.readline().split()
    fid.close()
    
    
    
        
    
    # Make Quals
    Quals = np.zeros((len(p), numSamples), dtype='uint') # initialize
    for i in range (numSamples):
        print('Loading quals matrix for sample: ' + str(i)) 
        print('Filename: ' + paths_to_quals_files[i]) 
        infile=sio.loadmat(paths_to_quals_files[i]) # from previous step, should include variable called p
        quals=infile['quals'].flatten()
        Quals[:,i]=quals[p]
    
    
    
    ## counts: counts for each base from forward and reverse reads at each candidate position for all samples
    print('Gathering counts data at each candidate position...\n')
    
    # Import list of directories for where to diversity.mat for each sample
    fname =  pwd + '/' + path_to_list_of_diversity_files  
    fid = open( fname, "r" ) # Input is space separated text file, in one line
    paths_to_diversity_files = fid.readline().split()
    fid.close()
    
    tempfile=sio.loadmat(paths_to_diversity_files[1]) 
    data=tempfile['data']
    size=np.shape(data)
    GenomeLength=size[1]
        
    # Make counts and coverage at the same time
    counts = np.zeros((8, len(p), numSamples),dtype='uint') # initialize
    all_coverage_per_bp = np.zeros((GenomeLength, numSamples),dtype='uint') # Added 2019.12.12
    indel_counter=np.zeros((2, len(p), numSamples), dtype='uint') # Added 2019.02.12
    for i in range (numSamples):
        print('Loading counts matrix for sample: ' + str(i)) 
        print('Filename: '+ paths_to_diversity_files[i]) 
        infile=sio.loadmat(paths_to_diversity_files[i]) 
        data=infile['data']
        counts[:,:,i]=data[0:8,p]
        all_coverage_per_bp[:,i]=sum(data[0:8,:])
        indel_counter[:,:,i]=data[38:40,p] # Added 2019.02.12 reads supporting indels and reads supporting deletions
    
    
    #print('Getting all the coverage information...\n')
    #[all_coverage_per_bp, ~, all_maf_per_bp] = get_all_coverage(SampleInfo, GenomeLength)
    
    # Normalize coverage by sample and then position; ignore /0 ; turn resulting inf to 0
    with np.errstate(divide='ignore',invalid='ignore'):
        array_cov_norm = ( all_coverage_per_bp - np.mean(all_coverage_per_bp,axis=0,keepdims=True) ) / np.std(all_coverage_per_bp,axis=0,keepdims=True) # ,keepdims=True maintains 2D array (second dim == 1), necessary for braodcasting
        array_cov_norm[ ~np.isfinite(array_cov_norm) ] = 0
        
        # 2nd normalisation
        array_cov_norm = ( array_cov_norm - np.mean(array_cov_norm,axis=1,keepdims=True) ) / np.std(array_cov_norm,axis=1,keepdims=True) # ,keepdims=True maintains 2D array (second dim == 1), necessary for braodcasting
        array_cov_norm[ ~np.isfinite(array_cov_norm) ] = 0
    

    ## turn into sparse csc matrices for more efficient computation
    # scale norm matrix by 1000 and save as int64 to slim matrix as much as possible
    all_coverage_per_bp_csc = sparse.csc_matrix(all_coverage_per_bp)
    array_cov_norm_scaled_csc = sparse.csc_matrix((np.round(array_cov_norm,3)*1000).astype('int64'))

    ## Save everything!
    print('Saving everything compressed...')
    
    with gzip.open(path_to_candidate_mutation_table, 'wb') as f: 
        pickle.dump([SampleNames, p, counts, Quals, in_outgroup, indel_counter], f)
    
    sparse.save_npz(path_cov_raw_sparse_matrix, all_coverage_per_bp_csc,compressed=True)
    sparse.save_npz(path_cov_norm_sparse_scale_matrix, array_cov_norm_scaled_csc,compressed=True)

    
    print('DONE')


if __name__ == "__main__":
    if len(sys.argv[1:])<8:
        print('Usage: Reads in 5 files and requires 3 output file names.')
    else:
        print('Grabbing args...')
        path_to_p_file=sys.argv[1]
        path_to_sample_names_file=sys.argv[2]
        path_to_outgroup_boolean_file=sys.argv[3]
        path_to_list_of_quals_files=sys.argv[4]
        path_to_list_of_diversity_files=sys.argv[5]
        path_to_candidate_mutation_table=sys.argv[6]
        path_cov_raw_sparse_matrix=sys.argv[7]
        path_cov_norm_sparse_scale_matrix=sys.argv[8]
        print('Calling main...')
        main(path_to_p_file, path_to_sample_names_file, path_to_outgroup_boolean_file, path_to_list_of_quals_files, path_to_list_of_diversity_files, path_to_candidate_mutation_table, path_cov_raw_sparse_matrix,path_cov_norm_sparse_scale_matrix)

