#########################################
# LIEBERMAN LAB SNAKEFILE FOR CASE STEP #
#########################################

# Version History:
# # 2019.02.12, Arolyn: rule candidate_mutation_table now calling make_candidate_mutation_table_snakemake_indels.m instead of make_candidate_mutation_table_snakemake.m
# # 2018.12, Felix/Arolyn: Original Snakefile from the lab hackathon

# Reminders:
# # Put your own email address in cluster.slurm.json so that you get emails about failures. No one else wants your slurm emails.


## USER defined variables (in theory do not need to be touched)
spls = "samples_case.csv"
SCRIPTS_DIRECTORY = "/scratch/mit_lieberman/scripts"
maxFQ = -30 # purity threshold (from mapping quality) for including position

## pre-snakemake
import sys, subprocess
sys.path.insert(0, SCRIPTS_DIRECTORY)

from read_samplesCSV import read_samples_caseCSV


''' PRE-SNAKEMAKE '''
## define couple of lists from samples_case.csv
## FORMAT: Path,Sample,ReferenceGenome,Outgroup
# NOTE: case analysis expects:
#	- only a single ref genome used for all samples
#	- unique sample IDs
#	- Col1==Paths, points to the snakemake directory used for mapping, not raw data! > only single entries accepted!
[PATH_ls, SAMPLE_ls, REF_Genome_ls, OUTGROUP_ls] = read_samples_caseCSV(spls)


''' SNAKEMAKE '''
rule all:
    input:
        # # Data links only # #
        # expand("data/vcf/{sampleID}_ref_{references}_outgroup{outgroup}.vcf.gz",zip,sampleID=SAMPLE_ls, references=REF_Genome_ls, outgroup=OUTGROUP_ls),
        # expand("data/qual/{sampleID}_ref_{references}_outgroup{outgroup}.quals.mat",zip,sampleID=SAMPLE_ls, references=REF_Genome_ls, outgroup=OUTGROUP_ls),
        # expand("data/diversity/{sampleID}_ref_{references}_outgroup{outgroup}.diversity.mat",zip,sampleID=SAMPLE_ls, references=REF_Genome_ls, outgroup=OUTGROUP_ls),
        # # Everything # #
        "2-candidate_mutation_table/candidate_mutation_table.mat",
        # # Including cleanup step # #
        # "logs/DONE_cleanUp",


rule build_data_links:
	output:
		vcf_links = expand("data/vcf/{sampleID}_ref_{references}_outgroup{outgroup}.vcf.gz",zip,sampleID=SAMPLE_ls, references=REF_Genome_ls, outgroup=OUTGROUP_ls),
		qual_mat_links = expand("data/qual/{sampleID}_ref_{references}_outgroup{outgroup}.quals.mat",zip,sampleID=SAMPLE_ls, references=REF_Genome_ls, outgroup=OUTGROUP_ls),
		div_mat_links = expand("data/diversity/{sampleID}_ref_{references}_outgroup{outgroup}.diversity.mat",zip,sampleID=SAMPLE_ls, references=REF_Genome_ls, outgroup=OUTGROUP_ls),
	log:
		"logs/build_links.log"
	run:
		import subprocess
		subprocess.run( "rm -fr data/ " ,shell=True) # clean it up prior run
		subprocess.run( "mkdir -p data/vcf/ data/qual/ data/diversity/ " ,shell=True)
		for i in range(len(SAMPLE_ls)):
			subprocess.run( "ln -fs -T " + PATH_ls[i] + "/6-diversity/" + SAMPLE_ls[i] + "_ref_" + REF_Genome_ls[i] + "*diversity.mat data/diversity/" + SAMPLE_ls[i] + "_ref_" + REF_Genome_ls[i] + "_outgroup" + OUTGROUP_ls[i] + ".diversity.mat" ,shell=True)
			subprocess.run( "ln -fs -T " + PATH_ls[i] + "/5-quals/" + SAMPLE_ls[i] + "_ref_" + REF_Genome_ls[i] + "*quals.mat data/qual/" + SAMPLE_ls[i] + "_ref_" + REF_Genome_ls[i] + "_outgroup" + OUTGROUP_ls[i] + ".quals.mat" ,shell=True)
			subprocess.run( "ln -fs -T " + PATH_ls[i] + "/4-vcf/" + SAMPLE_ls[i] + "_ref_" + REF_Genome_ls[i] + "*variant.vcf.gz data/vcf/" + SAMPLE_ls[i] + "_ref_" + REF_Genome_ls[i] + "_outgroup" + OUTGROUP_ls[i] + ".vcf.gz " ,shell=True)
			

rule variants2positions:
	input:
		variants = "data/vcf/{sampleID}_ref_{references}_outgroup{outgroup}.vcf.gz",
	params:
		REF_GENOME_DIRECTORY = "/scratch/mit_lieberman/reference_genomes/{references}/",
		outgroup_tag = "{outgroup}", # boolean (0==ingroup or 1==outgroup)
	log:
		"logs/variants2positions_{sampleID}_ref_{references}_outgroup{outgroup}.log"
	output:
		mat_positions = "1-temp_pos/{sampleID}_ref_{references}_outgroup{outgroup}_positions.mat",
	shell:
		"""
		module add mit/matlab/2015b; 
		matlab -r "path('{SCRIPTS_DIRECTORY}',path); generate_positions_single_sample_snakemake( '{input.variants}', '{output.mat_positions}', {maxFQ}, '{params.REF_GENOME_DIRECTORY}',{params.outgroup_tag}  ) " > {log} 2>&1		
		"""


rule combine_positions:
    input:
        mat_positions = expand("1-temp_pos/{sampleID}_ref_{references}_outgroup{outgroup}_positions.mat",zip,sampleID=SAMPLE_ls, references=REF_Genome_ls, outgroup=OUTGROUP_ls)
    params:
        file_other_p_to_consider = "add_positions/other_positions.mat",
        # REF_GENOME_DIRECTORY = "/scratch/mit_lieberman/reference_genomes/{references}/",
        REF_GENOME_DIRECTORY = expand("/scratch/mit_lieberman/reference_genomes/{references}/",references=set(REF_Genome_ls)), # expands to single reference genome!
    log:
    	"logs/combine_positions.log",
    output:
        string_input_p_positions = "1-temp_pos/string_file_other_p_to_consider.txt",
        mat_positions = "1-temp_pos/allpositions.mat",
    shell:
        """
        module add mit/matlab/2015b; 
        echo {input.mat_positions} > {output.string_input_p_positions}
        matlab -r "path('{SCRIPTS_DIRECTORY}',path); combine_positions_snakemake( '{output.string_input_p_positions}', '{params.file_other_p_to_consider}', '{output.mat_positions}', '{params.REF_GENOME_DIRECTORY}', {maxFQ} )" > {log} 2>&1	
        """


# build input for candidate_mutation_table
rule string_diversity_mat:
    input:
    	diversity_mat = expand("data/diversity/{sampleID}_ref_{references}_outgroup{outgroup}.diversity.mat",zip,sampleID=SAMPLE_ls, references=REF_Genome_ls, outgroup=OUTGROUP_ls),
    output:
    	string_diversity_mat = "1-temp_pos/string_diversity_mat.txt",
    run:
    	with open( output.string_diversity_mat ,"w") as f: 
    		print(*input.diversity_mat, sep=" ", file=f)


# build input for candidate_mutation_table
rule string_quals_mat:
    input:
    	quals_mat = expand("data/qual/{sampleID}_ref_{references}_outgroup{outgroup}.quals.mat",zip,sampleID=SAMPLE_ls, references=REF_Genome_ls, outgroup=OUTGROUP_ls),
    output:
    	string_qual_mat = "1-temp_pos/string_qual_mat.txt",
    shell:
    	"""
    	for i in {input.quals_mat}; do
    	 echo $i
    	done |tr '\n' ' ' |sed "s/ $//" > {output.string_qual_mat};
        """


# build input for candidate_mutation_table
rule string_sampleID_names:
    params:
    	sampleID_names = expand( "{sampleID}" , sampleID=SAMPLE_ls ),
    output:
    	string_sampleID_names = "1-temp_pos/string_sampleID_names.txt",
    shell:
    	"""
    	for i in {params.sampleID_names}; do
    	 echo $i
    	done |tr '\n' ' ' |sed "s/ $//" > {output.string_sampleID_names};
        """


# build input for candidate_mutation_table
rule string_outgroup_bool:
    params:
    	outgroup_bool = expand( "{outgroup}" , outgroup=OUTGROUP_ls ),
    output:
    	string_outgroup_bool = "1-temp_pos/string_outgroup_bool.txt",
    shell:
    	"""
    	for i in {params.outgroup_bool}; do
    	 echo $i
    	done |tr '\n' ' ' |sed "s/ $//" > {output.string_outgroup_bool};
        """


rule candidate_mutation_table:
    input:
    	mat_positions = "1-temp_pos/allpositions.mat",
        string_diversity_mat = "1-temp_pos/string_diversity_mat.txt",
    	string_qual_mat = "1-temp_pos/string_qual_mat.txt",
    	string_sampleID_names = "1-temp_pos/string_sampleID_names.txt",
    	string_outgroup_bool = "1-temp_pos/string_outgroup_bool.txt",
    output:
    	candidate_mutation_table = "2-candidate_mutation_table/candidate_mutation_table.npz",
    	cov_mat = "2-candidate_mutation_table/cov_mat.npz",
	norm_cov_mat = "2-candidate_mutation_table/norm_cov_mat.npz",
    log:
    	"logs/candidate_mutation_table.log",
    conda:
    	"envs/snakemake.yaml",
    shell:
    	"""
    	python3 build_candidate_mutation_table.py {input.mat_positions} {input.string_sampleID_names} {input.string_outgroup_bool} {input.string_qual_mat} {input.string_diversity_mat} {output.candidate_mutation_table} {output.cov_mat} {output.norm_cov_mat}
    	"""




rule cleanUp:
	input:
		candidate_mutation_table = "2-candidate_mutation_table/candidate_mutation_table.mat",
	params:
		temp_folder = "1-temp_pos"
	output:
		"logs/DONE_cleanUp"
	shell:
		" rm -rf {params.temp_folder} ; touch {output} "

