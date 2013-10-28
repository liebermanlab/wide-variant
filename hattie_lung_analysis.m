% Lung Specific Analysis

lung_analysis = 1; 

if lung_analysis == 1
    % plot pairwise distance within site 
    lung_pairwise_distances(mut_freq, SampleNames); 

    % plot number of unique isolates within and between sites 
    lung_get_all_clonality(mut_freq, SampleNames); 

    % plot distribution of pairwise distances
    [within, between, distance_matrix] = plot_all_isolate_pairwise_distances(mut_freq, SampleNames); 
end