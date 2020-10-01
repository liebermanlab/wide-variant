# -*- coding: utf-8 -*-
"""
Placeholder script
"""

print('Hello!')

import numpy as np
import pickle
import scipy.io as sio
import os
import sys
import gzip
from scipy import sparse

print('Done importing stuff...')

def main(path_to_p_file, path_to_sample_names_file, path_to_outgroup_boolean_file, path_to_list_of_quals_files, path_to_list_of_diversity_files, path_to_candidate_mutation_table, path_cov_raw_sparse_matrix,path_cov_norm_sparse_scale_matrix):

    print('Nothing to do...')


if __name__ == "__main__":
    if len(sys.argv[1:])<8:
        print('Usage: Reads in 5 files and requires 3 output file names.')
    else:
        path_to_p_file=sys.argv[1]
        path_to_sample_names_file=sys.argv[2]
        path_to_outgroup_boolean_file=sys.argv[3]
        path_to_list_of_quals_files=sys.argv[4]
        path_to_list_of_diversity_files=sys.argv[5]
        path_to_candidate_mutation_table=sys.argv[6]
        path_cov_raw_sparse_matrix=sys.argv[7]
        path_cov_norm_sparse_scale_matrix=sys.argv[8]
        main(path_to_p_file, path_to_sample_names_file, path_to_outgroup_boolean_file, path_to_list_of_quals_files, path_to_list_of_diversity_files, path_to_candidate_mutation_table, path_cov_raw_sparse_matrix,path_cov_norm_sparse_scale_matrix)

