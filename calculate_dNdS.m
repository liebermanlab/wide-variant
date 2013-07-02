function [ci_upper, ci_lower, dnds] = calculate_dNdS(annotation_mutgenes, cds, GenomeLength, ChrStarts, sequences)
    
    % Get observed number of N muts
    observed_N = div_count_observed_mut_types(annotation_mutgenes, 'N'); 
    observed_N_6types = div_matrix2_6types(observed_N); 
    
    % Get observed number of S muts
    observed_S = div_count_observed_mut_types(annotation_mutgenes, 'S'); 
    observed_S_6types = div_matrix2_6types(observed_S); 
    
    % Get null model
    if exist('steno_edited_null_dNdS.mat') == 2
        disp('Loading existing null probability model');
        null = load('steno_edited_null_dNdS'); 
        m = null.m; 
        m_coding_strand = null.m_coding_strand; 
        probN = null.probN; 
    else
        disp('Calculating null probability model... get a snickers'); 
        [m, m_coding_strand, probN] = div_mutation_type_probability_matrix(cds, GenomeLength, ChrStarts, sequences);
    end
    
    % Calculate expected number of N muts
    expected_N = dot(probN, observed_N_6types); 
    
    % Calculate 95% Binomial CI and dNdS ratio 
    [ci_upper, ci_lower, dnds] = binomialCIdNdS(observed_N_6types, observed_S_6types, expected_N); 
    
end