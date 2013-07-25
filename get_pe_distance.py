'''
HC 7/11/2013
Takes SAM file and saves all PE distsances 
USAGE: python get_pe_distance.py [sam file] 
''' 

import os, sys, math 

def parse_sam_file(fn): 
    # open file
    fh = open(fn, 'rU')
    
    distances = [] 
    counter = 0
    # parse SAM file
    for line in fh: 
        # ignore headers
        if line.startswith('@') is False: 
            lineout = line.split() 
            start_pos, end_pos = lineout[3], lineout[7] 
            if (start_pos is not 0) and (end_pos is not 0): 
                dist = math.fabs(int(end_pos) - int(start_pos)) 
                distances.append(dist) 
                counter += 1

    print 'Stored %i distances' %counter 

    # save distances matrix 
    output = open('pe_distances.csv', 'w') 
    output.write( ','.join(str(x) for x in distances) ) 
    output.close() 

if __name__ == '__main__': 
    filename = sys.argv[1]
    parse_sam_file(filename) 
