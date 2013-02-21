%% Set important variables each run

run_postfix='13_01_08'; %must match postfix in build_mutation_table_master.m

%run options
loadwindows=1;
promotersize=150;
ancestorsample=1;
referenceisancestor=0;

%filter parameters
qual_0=282;
strict_parameters=struct('minorfreqthreshold',.03, 'maxreads_perstrand_percentile',99,...
    'minreads_perstrand',30, 'minreads_perstrand_per_allele',2,'min_bq',19,'min_mq', 33, 'min_td', 7,...
    'max_td',93, 'max_sbp', 3,'max_percent_indels', .20, 'min_control_MAF', .98, ...
    'max_bqp', 200,'max_tdp',200, 'max_percent_ends', 1);


% ___ create parameters log ___ %

% add promoter size, qual_0 to log

log_parameters = strict_parameters; 
log_parameters.('promotersize') = promotersize; 
log_parameters.('qual_o') = qual_0; 

logfolder = 'log_analysis/'; 
if exist(logfolder, 'dir') ~= 7
    mkdir(logfolder)
end

save_structure_parameters(logfolder, log_parameters); 

% _____________________________ %



%% Set useful useful variables


path('../scripts',path)

NTs='ATCG';
scrsz = get(0,'ScreenSize');



%% Load output from build_mutation_table_master

SampleInfo = read_sample_names ;
SampleNames={SampleInfo(:).Sample}';
load(['mutation_table_' run_postfix])
load(['MutGenVCF_' run_postfix])
Nsample=numel(SampleNames);


cwindows=[]; fwindows=[];
if loadwindows==1
    load(['windows_' run_postfix])
end

showlegends=1;
if numel(SampleNames)>10
    showlegends=0;
end

%% Temporary:: Fix annotation structure
[geneloc,  cds, mutations, sequences] = annotate_mutations_auto_gb(positions,ScafNames,RefGenome) ;



%% Fixed mutations:  Quality threshold anaysis (set value for qual_0)

%Number of mutations
MutQual = ana_mutation_quality(Calls,Quals,1) ;
plot(qual_0,sum(MutQual>=qual_0),'dr', 'MarkerFaceColor', 'r', 'MarkerSize', 10)

%Pairwise distance between strains
step=floor(qual_0/10);
qual_th = sort([0:step:(max(MutQual)) qual_0 max(MutQual)-1]) ; 
figure(7);clf
set(gcf,'position',[10   570   560   420])
for m=1:length(qual_th)
    Calls_temp = Calls ; Calls_temp(Quals<qual_th(m)) = 'N' ;
    df = ana_strain_distance(Calls_temp(:,1:Nsample),SampleNames) ;
    title(sprintf('Quality > %g',qual_th(m))) ;
    %pause
end



%% Create useful matrices and vectors

%Ancestral nucleotide at each position
if ancestorsample>0
    [ancnt, ancnti] = ancestorfromsample(Calls,ancestorsample);
end
ancnt_m=repmat(ancnt,1,Nsample);
ancnti_m=repmat(ancnti,1,Nsample);

[maf, maNT, minorNT] = div_major_allele_freq(counts);

%Hasmutation
diversemutation=div_test_thresholds(counts,strict_parameters, coveragethresholds);
fixedmutation=((Calls~=ancnt_m) & (ismember(Calls,'ATCGDI')) & repmat(MutQual,1,Nsample)>=qual_0);
hasmutation= fixedmutation | diversemutation; %has mutation if diverse or call(from vcf file) ~= anct
minormutation=(hasmutation & (ancnti_m==maNT));

%Mutant allele frequency
[mutAF, mutantNT]=mutant_frequency(counts, hasmutation, ancnti, Calls);


%% Generate table

[q, annotation_all] = div_clickable_table(mutations, p, ancnti, counts,  fwindows, cwindows, mutAF, diversemutation, RefGenome, ScafNames, SampleInfo, ChrStarts, promotersize, showlegends);


%% More useful information

types=[annotation_all.type];
typesmatrix=repmat(types,1,Nsample);

genes=locustags2numbers({annotation_all.locustag});

    
%%

save(['mutation_analysis_' run_postfix], 'mutations', 'Nsample', 'mutAF', 'genes')






stop





%% Find coverable positions

if analyze_diversity==1
    fprintf(1,'Finding genomic positions with potential to call diversity\n');
    callablePos= div_genomic_positions_with_potential_to_call_diversity(SampleDirs, SampleNames, strict_parameters, coveragethresholds, GenomeLength, Parallel, jobsubmitoptions, generatehists);
    %remove from callable if diverse in control
    [maf, maNT, minorNT] = div_major_allele_freq(counts);
    callablePos(p(maf(:,1)<strict_parameters.min_control_MAF),:)=0;
end   

%coverage positions in isolates previously generated with a python script and stored
%in a text file (based on FQ in strain.vcf)

%need to write a better way to do this
