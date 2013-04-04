----- Overview -------------------------------------
This pipeline is intended to be run in 4 steps.
(A) Set up proper folders and documents.
(B) Process raw data files into sample-specific files, using 1 set of parameters per sample, or many (This is the 'experiment' step, build_experiment directories)
(C)  Gather raw data from every sample at every genomic position that was called as a potential variant in at least one sample  (This is part 1 of the 'case' step, build_mutation_table_master)
(D) Filter and produce interactive table and useful data structures (This is part 2 of the 'case' step, analysis_master)

The purpose of separating the experiment and case steps is, for each sample, to try many sets of filter and alignment parameters. Then we can select the best parameter set to move forward with when comparing across samples.


----- (A) Setup ----------------------------------------
1. Set up reference genome folder in Reference_Genomes. Add fasta and genebank files. Name the fasta file genome.fasta. There may be multiple genebank files-- name each genebank file the name of the corresponding chromosome/scaffold in the fasta file (e.g. NC_002695.1.gb).
2. Set up an experiment folder with samples.csv, read_filter_params.csv, alignment_params.csv. This is best in "/scratch" folder in orchestra AND CANNOT HAVE A SPACE IN THE PATH. See other folder for examples. In samples.csv, use spaces to indicate multiple alignment parameters for a given sample. Line delimiter must be written by excel (I think this is '\r' . Textwrangler doesn't seem to work).
3. Set up a case folder with sample_names.csv . If calling polymorphic positions, the isogenic control must be the first sample in sample_names.csv (Case folder). The last column of this folder may need to be adjusted after running part B. 


----- (B) Process raw data files ------------------
This step creates sam, pileup, vcf, and diversity.mat (info on this structure below) files for finding polymorphisms. It also creates a .bai file which helps for IGV viewer and samtools view, two ways of visualizing alignments. A future version of this pipeline should also create coverage maps at this point, as coverage maps do not depend on other samples. In short, any sample-specific file should be created in this step.
1. Open up MATLAB on interactive session on cluster in this experiment folder and run build_experiment_directories


----- (C) Gather raw data at potentially variant sites---------
This step creates 3 structures which summarize the position at every potentially variant site across all samples. Descriptions below.
1. Set up a case folder with sample_names.csv . If calling polymorphic positions, the isogenic control must be the first sample in sample_names.csv (Case folder). 
2. When setting up case folder, each sample must have a single alignment parameter chosen. Each case can only have one Reference Genome but can have different parameter sets. You can chose the parameter sets by using the matlab function choose_alignment_parameters([path]/alignment_summary)
3. Inspect header of build_mutation_table_master for hardcoded variables (logged each run)
4. Open up MATLAB on interactive session on cluster in this case folder and run build_mutation_table_master

--- (D) Filter and produce interactive table and useful data structures -----
1. Inspect header of analysis_master for hardcoded variables (logged each run). Ensure value of postfix matches what was run for step ©.
2. Open up MATLAB locally in this the folder and run analysis_master





-----Example---------------------------------------
(Make reference genome folder with fasta and genebank files)
ssh USER@orchestra.med.harvard.edu          (on HMS machine, simply: ssh orchestra)
(Make experiment folder and csv files in scratch)
cd [experimentfolder]
bsub -Is -q sysbio_int_2d matlab
scriptsdirectory =[scriptsdirectory]
path(scriptsdirectory,path)
build_experiment_directories(scriptsdirectory)


(Make case folder and csv file)
cd [casefolder]
build_mutation_table_master


(switch to running matlab locally, from Case Folder)

analysis_master



---About mutation_table_[postfix].mat--------------
This contains a lot of other variables.

The following variables are straightforward
:::RefGenome-- name of reference 
:::ScafNames -- names of chromosomes/scaffolds of reference genome
:::GenomeLength -- total length of reference genome
:::cds -- A structure containing some basic information about each coding sequence on the reference genome (this is generally not used)
:::sequences -- The actual nucleotide sequence of the reference

The following variables require are not as simple:::
:::positions -- List of genomic positions with potential variants. Positions is size [2] x [npositions], with the first column being the chromosome/scaffold number, and the second column the position on that scaf
:::p -- A vectorized index of positions. p and positions can be turned from one to the other using p2chrpos(p,ChrStarts) and chrpos2index(positions,ChrStarts) 
:::ChrStarts -- A list of where in the vector 1:GenomeLength each chromosome/scaffold starts. 

:::geneloc -- The cds number corresponding to each variant position.
:::mutations-- A structure which contains information about the annotation of each genomic position irrespective of what mutations were actually called there. It includes the gene mutated, if applicable, the reference nt and AA sequences of the gene carrying the mutation. It also includes the very useful mutations.AA, which indicates what amino acid is encoded at this codon if there is an A, T, C or G, in this position (e.g. 'mutations(3).AA='FIVL'). It assumes no other mutations relative to the reference in this codon. 

:::Calls -- The best guess of a NT call at this each position for each strain. Size [npositions] x [nstrains].
:::Quals -- The FQ score (consensus) produced by samtools at each position for each strain. Size [npositions] x [nstrains].

:::Counts -- The same 39 attributes in diversity.mat, but instead of all genomic positions for a single strain, it is at all variant positions in all strains. Size 39 x [npositions] x [nstrains].




---About windows_[postfix].mat--------------
This is a matrix contains cwindows and fwindows.  windows stores coverage information and fwindows stores major allele frequency information. Each matrix holds the coverage or frequency of the genomic positions in the +/- window_size (a parameter designated at the top of build_mutation_table_master. Generally set to 500bp, producing a 1kb window) around each variant position for each sample. Each is a 3d matrix of size [(window_size*2)+1] x [npositions] x [nsamples].

---About MutGenVCF_[postfix].mat--------------
This data structure may actually not need to be passed on to analysis_master. This is a data structure which stores all of the information in the samtools-produced .vcf files, for all samples, at each potentially variant position. The FQ scores ('consensus scores') are generally what is used in downstream analysis.



---About diversity.mat--------------
This is a 39 x [Genomelength] structure that contains, for each position on the genome
1-4: Number of reads aligning to the Forward strand supporting A,T,C, and G
5-8: Number of reads aligning to the Reverse strand supporting a,t,c, and g
9-16: Average phred base qualities of calls supporting A,T,C, and G
17-24: Average mapping quality of reads supporting A,T,C, and G
24-32: Average tail distance of reads supporting A,T,C, and G
33: -log10(Fishers exact test - probability of counts (number of reads) supporting major and minor alleles in forward and reverse direction coming from the same distribution)
34: -log10(Ttest - probability of base qualities of major and minor alleles coming from the same distribution)
35: -log10(Ttest - probability of mapping qualities of major and minor alleles coming from the same distribution)
36: -log10(Ttest - probability of tail distances on the forward strand for major and minor alleles coming from the same distribution)
37: -log10(Ttest - probability of  tail distances on the reverse strand of major and minor alleles coming from the same distribution)
38: number of calls aligning that are at either end of a red
39: number of reads supporting an within 3 bp of this position
