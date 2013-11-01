%% Variables you might want to change 

%path('/Users/tdl7/Dropbox/illumina_pipeline',path)



global ANALYZE_DIVERSITY; global qual_0; global CONTROLSAMPLE; global LABEL_SIZE;
SampleInfo = read_sample_names ;



%important options
qual_0=130; % want to use FQ and not qual from VCF 
goodsamples=1:numel(SampleInfo);

%run options
run_postfix='13_04_02'; %must match postfix in build_mutation_table_master.m

%diversity
ANALYZE_DIVERSITY= 0;
CONTROLSAMPLE=1; % deep isogenic control 

%other options
onlySNPs=1; 
promotersize=150;

%choose ancestor setting
referenceisancestor=1;
ancestoriscontrol=0;
ancestorismode=0; %use mode of major alleles as ancestor


%Display options
LABEL_SIZE=2;
loadwindows=1; % f and c windows 
barcharttype=3; %1 means shows all strains; 2 means show 2 used for MutQual, 3 means is same #2, but also one that was clicked




%% Initialize
analysis_master_initialize; 




make_useful_matrices; 







%% Generate table -- inspect lower MutQuals and toggle qual_0

QualSort=1;
QualCutOff=1;
[annotation_all, sorted_table_data] = div_clickable_table_isolate_calls(mutations, Calls, p, ancnti, ...
                                            counts,  fwindows, cwindows, ...
                                            hasmutation, MutQual, MutQualIsolates, ...
                                            RefGenome, ScafNames, SampleInfo, ...
                                            ChrStarts, promotersize, showlegends, ...
                                            QualSort, QualCutOff, qual_0);   



%% Clustergram

% OPTION 1. 
% filled_annotation = fill_in_annotation(annotation_all); 

% OPTION 2. 
% for specifying whether NS, Syn, P, I mutation in clustergram 
% filled_annotation = cell(numel(annotation_all),1);
% for i = 1:numel(filled_annotation)
%     filled_annotation{i} = [annotation_all(i).type ' | ' annotation_all(i).annotation];
% end

% MAKE CLUSTERGRAM.
% co = clustergram(mut_freq);
% set(co,'ColumnLabels', SampleNames, 'RowLabels', filled_annotation, ...
%         'Colormap', redbluecmap);

%% dNdS

%[ci_u, ci_l, dnds] = calculate_dNdS(annotation_all, cds, GenomeLength, ChrStarts, sequences); 


%% Find genes not unaligned to

% get_missing_genes; 

%% cooccurrence = plot_isolates_covariance(mut_freq, annotation_mutgenes); 

%% Get number of genes mutated more than once

% get_multiple_mutated_genes; 

%% Generate phylogeny 

make_phylogeny; 
 
%% Save

% save(['mutation_analysis_' run_postfix], 'mutations', 'Nsample', 'mutAF', 'genes')
