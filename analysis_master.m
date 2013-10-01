% Initialize 
analysis_master_initialize; 
make_useful_matrices; 

%% Generate table -- inspect lower MutQuals and toggle qual_0

QualSort=0;
QualCutOff=1;
[annotation_all, sorted_table_data] = div_clickable_table_isolate_calls(mutations, Calls, p, ancnti, ...
                                            counts,  fwindows, cwindows, ...
                                            hasmutation, MutQual, ...
                                            RefGenome, ScafNames, SampleInfo, ...
                                            ChrStarts, promotersize, showlegends, ...
                                            QualSort, QualCutOff, qual_0);                                        

%% Fill in missing annotation

% filled_annotation = fill_in_annotation(annotation_all); 

filled_annotation = cell(numel(annotation_all),1); 
for i = 1:numel(filled_annotation)
    filled_annotation{i} = [annotation_all(i).type ' | ' annotation_all(i).annotation];
end



%% Clustergram
mut_freq = get_fixed_mut_freq(Calls, MutQual, qual_0, ancnti);
co = clustergram(mut_freq);
set(co,'ColumnLabels', SampleNames, 'RowLabels', filled_annotation, ...
        'Colormap', redbluecmap);

%% Lung Specific Analysis

% plot pairwise distance within site 
lung_pairwise_distances(mut_freq, SampleNames); 

% plot number of unique isolates within and between sites 
lung_get_all_clonality(mut_freq, SampleNames); 

% plot distribution of pairwise distances
plot_all_isolate_pairwise_distances(mut_freq, SampleNames); 

%% dNdS

% [ci_u, ci_l, dnds] = calculate_dNdS(annotation_genes, cds, GenomeLength, ChrStarts, sequences); 


%% Find genes not unaligned to

% get_missing_genes; 

%% cooccurrence = plot_isolates_covariance(mut_freq, annotation_mutgenes); 

%% Get number of genes mutated more than once

% get_multiple_mutated_genes; 

%% Generate phylogeny 

% make_phylogeny; 
 
%% Save

% save(['mutation_analysis_' run_postfix], 'mutations', 'Nsample', 'mutAF', 'genes')
