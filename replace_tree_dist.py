import os, sys, re, math

# param
distscale = 0.00901

def replace_dist(line):
    pattern = '(?<=:)\d+.\d+'
    # calculate values to substitute with 
    matches = re.findall(pattern, line) 
    m_sub = [str(int(math.ceil(float(m)/distscale))) for m in matches] 
    
    # replace 
    subs = iter(m_sub) 
    newline = re.sub(pattern, lambda _: next(subs), line) 
    return newline

def main(in_file, out_file): 
    fh = open(in_file, 'rU')
    output = open(out_file, 'w') 
    for line in fh:
        # get new line to write
        newline = replace_dist(line)     
        output.write(newline) 
    output.close()

if __name__ == '__main__': 
    if len(sys.argv[1:]) is not 2: 
        print 'Please specify input tree and output file name' 
        sys.exit(1)
    inputfile = sys.argv[1] 
    outfile = sys.argv[2] 
    main(inputfile, outfile) 
