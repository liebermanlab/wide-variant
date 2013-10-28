function build_mutation_table_master(scriptsdirectory, CLUSTERDIR)

if nargin < 1 
     scriptsdirectory='/home/hc168/illumina_pipeline'; 
end

if nargin < 2
    CLUSTERDIR = '/groups/kishony';
end


%% Important variables to set each time

run_postfix='13_04_02';
TEMPFOLDER = '/hms/scratch1/hattie'; 
 

%Run in cluster?
Parallel=1;
onlysnps=1;
jobsubmitoptions1='sysbio_2h'; %short -W 0:15
jobsubmitoptions2='sysbio_12h'; %short -W 0:15

global RUN_ON_CLUSTER; RUN_ON_CLUSTER = 1;


positionfiles={}; %get data from other positions? (found in analysis of single isolates, previous study)

gscatters=0; %output options -- this useful graph has not been tested since updating code

window_size=500; %parameters used during initial data structure generation

analyze_diversity=1; %analyze diversity?

looseFQmax=-30;

loose_parameters=struct('minorfreqthreshold',.01, 'minreads_perstrand',2,...
    'maxreads_perstrand_percentile', 100,'minreads_perstrand_per_allele',1,...
    'min_bq',15,'min_mq' ,30, 'min_td', 10, 'max_td',90, 'max_sbp', 5,...
    'max_bqp', 255,'max_tdp',255, 'max_percent_ends', .50, 'max_percent_indels', .30, 'min_control_MAF', .01);

%most threshold checks are strictly > or strictly <
%loose parameters doesn't have an upper coverage thershold yet

% min_td and max_td are not symmetrical relative to read length of 100
% because some reads were trimmed prior to alignment
%maxreads_perstrand_percentile is which threshold in list .01:.01:1 ...
%    e.g. 98 is 98 percentile of covered positions



%% create parameters log 

% add window size to log 
log_parameters = loose_parameters; 
log_parameters.('window_size') = window_size; 

% set log directory
logfolder = 'log_build_mutation_table'; 
if exist(logfolder, 'dir') ~= 7
    mkdir(logfolder);
end

global TEMPORARYFOLDER;
timestamp=save_structure_parameters(logfolder, log_parameters);
TEMPORARYFOLDER=[TEMPFOLDER timestamp];

mkdir(TEMPORARYFOLDER)


%% Set path

global SCRIPTSPATH;

%For use as script
SCRIPTSPATH = scriptsdirectory;

%For use as function
% if nargin > 0
%      SCRIPTSPATH = scriptspath;
% else
%    a=pwd; SCRIPTSPATH=[a(1:find(double(a)==47,1,'last')) 'scripts']; 
%    if ~exist(SCRIPTSPATH, 'dir')
%        fprintf(['Could not find ' SCRIPTSPATH '...\n'])
%        error('Error: Must run from an experiment folder, with scripts folder in parent directory')
%    else
%        path(SCRIPTSPATH,path);
%    end
% end

fprintf(['Usings scripts directory: ' SCRIPTSPATH  '\n']);



%% Set main folder
if RUN_ON_CLUSTER == 1
    mainfolder=CLUSTERDIR;
else
    mainfolder='/Volumes/sysbio/KISHONY LAB/illumina_pipeline';
end



%% Read in csv files 


% Get SampleInfo and ScafNames
SampleInfo = read_sample_names ;

SampleNames={SampleInfo(:).Sample}';

NSample = length(SampleInfo) ;
RefGenome = {} ;
SampleDirs = {} ;

% Ensure that there is one reference genome
for i=1:NSample
    SampleDirs{i} = [SampleInfo(i).ExperimentFolder '/' SampleInfo(i).Sample '/' SampleInfo(i).AlignmentFolder ] ;
    ainfo = load([SampleDirs{i} '/alignment_info']) ;
    RefGenome{i} = ainfo.ae.Genome ;
end




RefGenome = unique(RefGenome) ;
if length(RefGenome)>1
    error('Must compare samples aligned to the same reference genome')
end

% Get Scafold names and build useful informations
RefGenome = RefGenome{1} ;
fr = fastaread([mainfolder '/Reference_Genomes/' RefGenome '/genome.fasta']) ;
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

cp = generate_positions(SampleDirs, SampleNames, GenomeLength, ScafNames, ChrStarts, looseFQmax, onlysnps, 100000, Parallel, jobsubmitoptions1);
fprintf(1,'Found %g positions where samtools called a variant in at least one sample \n',length(cp)) ;


op=[];
if numel(positionfiles)>0
    for i=1:numel(positionfiles)
        other=load(positionfiles{i});
        op=[op; chrpos2index(other.Positions, ChrStarts)];
    end
end

dp=[];
if analyze_diversity
    fprintf(1,'\nFinding single nucleotide positions with within-sample polymorphisms...');
    
    %Find diverse positions
    [dp, numfields, coveragethresholds] = find_diverse_positions(loose_parameters, SampleDirs, SampleNames, gscatters, Parallel, jobsubmitoptions1);

    fprintf(1,'Found %i positions with within-sample polymorphism that meets loose parameters in at least 1 sample \n',length(dp)) ;
    fprintf('\nNumber of samples %i\n', numel(SampleNames)); 
    memreq=2*4*length(dp)*numel(SampleNames)*(2*window_size)/(10^6);
    fprintf(1,['Please ensure that enough memory was requested when starting matlab session (use -R rusage[mem=' num2str(memreq) ']) (memory in MB) \n']);
    fprintf(1,'If p is large, frequency and coverage windows are not generated-- use smaller window or stricter parameters\n');
end 


allp=unique([dp; op; cp;]);
p=sort(allp);
p=p(p>0); %remove any 0s
positions=p2chrpos(p,ChrStarts);


%% Get counts and mutgenvcf

fprintf(1,'\nAcquiring detailed information at each potential position...\n');
fprintf(1,'vcf...');
[Calls, Quals] = gather_vcf_info(positions,SampleDirs,SampleNames, ScafNames,ChrStarts,jobsubmitoptions2,Parallel) ;                


fprintf(1,'diversity.mat...');
[counts, fwindows, cwindows] = generate_diversity_struct(SampleDirs, SampleNames, p, numfields, window_size, Parallel, jobsubmitoptions2) ;


%% Copy over select files so they don't get erased from scratch

fprintf(1,'\n\nCopying files over from scratch to pwd for storage...\n');

cmds={}; 
for i=1:NSample
    if ~exist(SampleNames{i},'dir')
        %don't use run_parallel_unix_commands_fast, because we want this to
        %run in background, and we aren't using the output at all
        mkdir(SampleNames{i})
    end
    cmds{end+1}=['cp ' SampleDirs{i} '/* ' SampleNames{i} '/'];
end
run_parallel_unix_commands_fast(cmds,jobsubmitoptions1,2);



%% Annotations

%This step generates an extra data structure containing information about
%the genomic position mutations -- can take > 10 minutes
[geneloc,  cds, mutations, sequences] = annotate_mutations_auto_gb(positions,ScafNames,RefGenome) ;


save(['mutation_table_' run_postfix], 'RefGenome', 'ScafNames', 'ChrStarts', 'GenomeLength', 'p', 'positions', 'coveragethresholds', 'counts',  'geneloc', 'cds', 'mutations', 'Calls', 'Quals', 'sequences', '-v7.3')
save(['windows_' run_postfix], 'fwindows', 'cwindows', '-v7.3')
% save(['MutGenVCF_' run_postfix], 'MutGenVCF', '-v7.3')

end

