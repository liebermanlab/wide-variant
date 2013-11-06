import os, sys, subprocess

# params
qname = 'sysbio_12h'; 
pathsdir = '/home/hc168/illumina_pipeline'
outputdir = 'unaligned_contig/'
hash_length = 31
min_contig_length = 5000

if __name__ == '__main__': 
    # Step 1. Merge all unaligned reads into 2 master files (FW + RV) 
    cmds = ["bsub", "-o", "combine_unaligned_out.txt", "-e", "combine_unaligned_err.txt", "-q", qname]
    cmds.append('"python /home/hc168/illumina_pipeline/build_contig_from_unaligned_reads.py"')
    subprocess.call(cmds)
    print '\nSubmitted job to combine unaligned reads' 
    # check that it finished 
    # search for words "finished" in combine_unaligned_out.txt 

    # Step 2. Run velveth
    velveth_cmds = [pathsdir+"/velvet_1.2.10/velveth", outputdir, str(hash_length), 
            "-fastq", "-shortPaired", "unaligned.all.1.fastq", "unaligned.all.2.fastq"]
    print '\nCalling velveth'
    subprocess.call(velveth_cmds) 
    
    # Step 3. Run velvetg 
    velvetg_cmds = [pathsdir+"/velvet_1.2.10/velvetg", outputdir, "-cov_cutoff", "10", "-min_contig_lgth", str(min_contig_length)]
