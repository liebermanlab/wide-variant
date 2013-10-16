%% Set important variables each run

global CONTROLSAMPLE


run_postfix='13_04_02'; %must match postfix in build_mutation_table_master.m

%run options
onlySNPs=1; 
loadwindows=0; % f and c windows 
promotersize=150;
CONTROLSAMPLE=1; % deep isogenic control 
referenceisancestor=0;
ancestoriscontrol=1;
ancestorismode=1; %use mode of major alleles as ancestor
qual_0=60; % want to use FQ and not qual from VCF 

MutQualFig=0; 
lung_analysis=0; 

strict_parameters = struct( 'minorfreqthreshold',           .03, ...
                            'maxreads_perstrand_percentile', 99, ...
                            'minreads_perstrand',            10, ...
                            'minreads_perstrand_per_allele', 4,...
                            'min_bq',                        19,...
                            'min_mq',                        33, ...
                            'min_td',                        7, ... % avg tail dist of each allele
                            'max_td',                        93, ...
                            'max_sbp',                       3,...  % p-val fisher's exact test of strand bias
                            'max_percent_indels',           .20, ...
                            'min_control_MAF',              .985, ...
                            'max_bqp',                      200, ...% t-test whether qual diff between two alleles
                            'max_tdp',                      200, ...% t-test whether tail dist diff b/t two alleles
                            'max_percent_ends',               1 ... 
                        );

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



%% Load output from build_mutation_table_master

SampleInfo = read_sample_names ;
SampleNames={SampleInfo(:).Sample}';
load(['mutation_table_' run_postfix])
% load(['MutGenVCF_' run_postfix])




%% Set useful useful variables


NTs='ATCG';
scrsz = get(0,'ScreenSize');
if onlySNPs
     acceptabletypes=NTs;
     Calls(~ismember(Calls,'ATCG'))='N';
else
     acceptabletypes='ATCGDI';
end


Nsample=numel(SampleNames);


cwindows=[]; fwindows=[];
if loadwindows==1
    load(['windows_' run_postfix])
end

showlegends=1;
if numel(SampleNames)>10
    showlegends=0;
end

%% Fixed mutations:  Quality threshold anaysis (set value for qual_0)

%Number of mutations
[MutQual, MutQualIsolates] = ana_mutation_quality(Calls,Quals,MutQualFig) ;
% plot(qual_0,sum(MutQual>=qual_0),'dr', 'MarkerFaceColor', 'r', 'MarkerSize', 10)

%Pairwise distance between strains
step=floor(qual_0/10);
qual_th = sort([0:step:min([max(MutQual)-1 qual_0+20]) qual_0]) ; 