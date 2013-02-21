%% Important variables to set each time

run_postfix='13_01_08';


%Run in cluster?
Parallel=1;
jobsubmitoptions='sysbio_2h'; %short -W 0:15

%get data from other positions? (found in analysis of single isolates, previous study)
get_data_at_other_specific_positions=0;
positionfiles={};

%output options -- this useful graph has not been tested since updating code
gscatters=0;

%parameters used during initial data structure generation
window_size=500;

%analyze diversity?
analyze_diversity=1;


loose_parameters=struct('minorfreqthreshold',.01, 'minreads_perstrand',10,...
    'maxreads_perstrand_percentile', 100,'minreads_perstrand_per_allele',1,...
    'min_bq',15,'min_mq' ,30, 'min_td', 10, 'max_td',90, 'max_sbp', 5,...
    'max_bqp', 255,'max_tdp',255, 'max_percent_ends', .50, 'max_percent_indels', .30, 'min_control_MAF', .01);

%most threshold checks are strictly > or strictly <
%loose parameters doesn't have an upper coverage thershold yet

% min_td and max_td are not symmetrical relative to read length of 100
% because some reads were trimmed prior to alignment
%maxreads_perstrand_percentile is which threshold in list .01:.01:1 ...
%    e.g. 98 is 98 percentile of covered positions

% ___ create parameters log ___ %

% add window size to log 
log_parameters = loose_parameters; 
log_parameters.('window_size') = window_size; 

% set log directory
logfolder = 'log_build_mutation_table'; 
if exist(logfolder, 'dir') ~= 7
    mkdir(logfolder)
end

save_structure_parameters(logfolder, log_parameters); 

% _____________________________ %



%% Read in csv files 
path('../scripts',path)


% Get SampleInfo and ScafNames
SampleInfo = read_sample_names ;

SampleNames={SampleInfo(:).Sample}';

NSample = length(SampleInfo) ;
RefGenome = {} ;
SampleDirs = {} ;

% Ensure that there is one reference genome
for i=1:NSample
    SampleDirs{i} = ['../' SampleInfo(i).ExperimentFolder '/' SampleInfo(i).Sample '/' SampleInfo(i).AlignmentFolder ] ;
    ainfo = load([SampleDirs{i} '/alignment_info']) ;
    RefGenome{i} = ainfo.ae.Genome ;
end




RefGenome = unique(RefGenome) ;
if length(RefGenome)>1
    error('Must compare samples aligned to the same reference genome')
end

% Get Scafold names and build useful informations
RefGenome = RefGenome{1} ;
fr = fastaread(['../Reference_Genomes/' RefGenome '/genome.fasta']) ;
GenomeLength=0;
ChrStarts=[];

ScafNames = {fr.Header} ;
for i=1:length(ScafNames)
    f=find(ScafNames{i}==' ',1) ;
    ScafNames{i} = ScafNames{i}(1:f-1) ;
    ChrStarts(end+1)=GenomeLength;
    GenomeLength=GenomeLength+numel(fr(i).Sequence);
end




%% Get all positions

fprintf(1,'\n\nFinding positions with at least 1 fixed mutation...\n');

cp = generate_positions(SampleDirs, SampleNames, ScafNames, ChrStarts, 100000, Parallel, jobsubmitoptions);
fprintf(1,'Found %g positions where samtools called a variant in at least one sample \n',length(cp)) ;


op=[];
if get_data_at_other_specific_positions==1
    
    for i=1:numel(positionfiles)
        other=load(positionfiles{i});
        op=[op; chrpos2index(other.Positions, ChrStarts)];
    end
end

dp=[];
if analyze_diversity
    fprintf(1,'\nFinding single nucleotide positions with diversity...');
    if ~exist('intermediate_diversity_matfiles')
        mkdir('intermediate_diversity_matfiles')
    end
    
    %Find diverse positions
    [dp, numfields, coveragethresholds] = find_diverse_positions(loose_parameters, SampleDirs, SampleNames, gscatters, Parallel, jobsubmitoptions);

    fprintf(1,'Found %g diverse positions that meet loose parameters in at least 1 sample \n',length(dp)) ;
    memreq=2*4*length(dp)*numel(SampleNames)*(2*window_size)/1000;
    fprintf(1,['Ensure that enough memory was requested when starting matlab session (use -R rusage[mem=' num2str(memreq) '])\n']);
    fprintf(1,'If p is large, frequency and coverage windows are not generated-- use smaller window or stricter parameters\n');
end 


allp=unique([dp; op; cp;]);
p=sort(allp);
positions=p2chrpos(p,ChrStarts);


%% Get counts and mutgenvcf

fprintf(1,'\nAcquiring detailed information at each potential position...\n');
fprintf(1,'vcf...');
[MutGenVCF, Calls, Quals] = generate_mutgenvcf_auto(positions,SampleDirs,SampleNames, ScafNames,jobsubmitoptions,Parallel) ;                
fprintf(1,'diversity.mat...');
[counts, fwindows, cwindows] = generate_diversity_struct(SampleDirs, SampleNames, p, numfields, window_size, Parallel, jobsubmitoptions) ;


%% Annotations

%This step generates an extra data structure containing information about
%the genomic position mutations -- can take > 10 minutes
[geneloc,  cds, mutations, sequences] = annotate_mutations_auto_gb(positions,ScafNames,RefGenome) ;


save(['mutation_table_' run_postfix], 'RefGenome', 'ScafNames', 'ChrStarts', 'GenomeLength', 'p', 'positions', 'coveragethresholds', 'counts',  'geneloc', 'cds', 'mutations', 'Calls', 'Quals', 'sequences', '-v7.3')
save(['windows_' run_postfix], 'fwindows', 'cwindows')
save(['MutGenVCF_' run_postfix], 'MutGenVCF')



