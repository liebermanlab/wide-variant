%% Set important variables each run

MutQualFig=1; 

STRICT_PARAMETERS = struct( 'minorfreqthreshold',           .10, ...
                            'maxreads_perstrand_percentile', 99, ...
                            'minreads_perstrand',            10, ...
                            'minreads_perstrand_per_allele', 2,...
                            'min_bq',                        19,...
                            'min_mq',                        33, ...
                            'min_td',                        7, ... % avg tail dist of each allele
                            'max_td',                        93, ...
                            'max_sbp',                       3,...  % p-val fisher's exact test of strand bias
                            'max_percent_indels',           .20, ...
                            'min_control_MAF',              .95, ...
                            'max_bqp',                      200, ...% t-test whether qual diff between two alleles
                            'max_tdp',                      200, ...% t-test whether tail dist diff b/t two alleles
                            'max_percent_ends',               1, ... 
                            'max_mqp',                      200 ...% t-test whether qual diff between two alleles
                        );

% ___ create parameters log ___ %

% add promoter size, qual_0 to log

log_parameters = STRICT_PARAMETERS; 
log_parameters.('promotersize') = promotersize; 
log_parameters.('qual_o') = qual_0; 

logfolder = 'log_analysis/'; 
if exist(logfolder, 'dir') ~= 7
    mkdir(logfolder)
end

save_structure_parameters(logfolder, log_parameters); 

% _____________________________ %




%% Set useful useful variables


NTs='ATCG';
scrsz = get(0,'ScreenSize');
if onlySNPs
     acceptabletypes=NTs;
     Calls(~ismember(Calls,'ATCG'))='N';
else
     acceptabletypes='ATCGDI';
end


cwindows=[]; fwindows=[];
if loadwindows==1
    load(['windows_' run_postfix])
end

showlegends=1;
if numel(SampleNames)>10
    showlegends=0;
end


SampleNames=SampleNames(goodsamples);
counts=counts(:,:,goodsamples);
cwindows=cwindows(:,:,goodsamples);
fwindows=fwindows(:,:,goodsamples);
Quals=Quals(:,goodsamples);
Calls=Calls(:,goodsamples);
%dindel_freqs=dindel_freqs(:,goodsamples);
%dindel_ids=dindel_ids(:,goodsamples);
%dindel_quals=dindel_quals(:,goodsamples);

Nsample=numel(SampleNames);


