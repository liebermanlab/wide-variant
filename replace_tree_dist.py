import os, sys, re, math

# param
lung_distscale = 0.00901

def replace_dist(line, distscale):
    pattern = '(?<=:)\d+.\d+'
    # calculate values to substitute with 
    matches = re.findall(pattern, line) 
    m_sub = [str(int(math.floor(float(m)/float(distscale)))) for m in matches]
    
    # replace 
    subs = iter(m_sub) 
    newline = re.sub(pattern, lambda _: next(subs), line) 
    return newline

def main(in_file, out_file, dist_scale): 
    fh = open(in_file, 'rU')
    output = open(out_file, 'w') 
    for line in fh:
        # get new line to write
        newline = replace_dist(line, dist_scale)     
        output.write(newline) 
    output.close()

if __name__ == '__main__': 
    '''
    USAGE: python replace_tree_dist.py [infile] [outfile] [distscale]
    '''
    if len(sys.argv[1:]) < 2: 
        print 'Please specify input tree and output file name' 
        sys.exit(1)
    elif len(sys.argv[1:]) is 2: 
        distscale = lung_distscale

    inputfile = sys.argv[1] 
    outfile = sys.argv[2] 
    distscale = sys.argv[3]
    main(inputfile, outfile, distscale) 
