function [genome_calls_all] = get_genome_counts_all(samplenames, genome_size) 
    % 2013 by Hattie Chung 

    % USAGE: get_all_covered_genes('SA*')
    
    % RETURNS: 
    % ___ all_coverage: [4 x GenomeLength] calls at each bp
    
    % params 
    threshold = 1; 
    genome_calls_all = uint16(zeros(4, length(samplenames), genome_size));
    
    % logical strucuture
    % 1 if position is found, 0 if missing
    for i = 1:length(samplenames)    
        sname = samplenames{i};
        % get covered positions
        fprintf('\nGetting coverage for isolate %i of %i\n', i, length(samplenames))
        [coverage] = get_genome_counts_sample(sname, threshold); 
        genome_calls_all(:,i,:) = coverage;
    end
end