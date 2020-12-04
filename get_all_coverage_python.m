function [all_coverage_per_bp, cov_modes, all_maf_per_bp] = get_all_coverage_python(sampleinfo, genome_size) 
    % NOTE: Alex: 10/18/20 this function is only for Dolosa samples and is meant to
    % quickly deal with the memory issues occuring with this project
    % 
    % OUTPUT: where n is number of samples 
    % all_coverage_per_bp: [n x GenomeLength] coverage at each bp
    % cov_modes: [n x 1] coverage mode for each sample  
    
    % params 
    all_coverage_per_bp = uint16(zeros(length(sampleinfo), genome_size));
    all_maf_per_bp = uint16(zeros(length(sampleinfo), genome_size));
    cov_modes = zeros(length(sampleinfo),1); 
    % logical strucuture
    % 1 if position is found, 0 if missing
    for i = 1:length(sampleinfo)   
        samplein = sampleinfo{i};
        filename = ['/scratch/mit_lieberman/projects/gga_bdelosa/mapping_python_10_09_20/6-diversity/', samplein, '_ref_Bdolosa_new_aligned.sorted.strain.variant.diversity.mat'];
        % get covered positions
        fprintf('\nGetting coverage for isolate %i of %i\n', i, length(sampleinfo))
        [coverage, coverage_mode, maf] = get_sample_coverage(filename, genome_size);  %fprintf(1,['\n' num2str(numel(coverage)) ' ' num2str(genome_size)]);
        all_coverage_per_bp(i,:) = coverage;
        all_maf_per_bp(i,:)=1000*maf;
        cov_modes(i) = coverage_mode(1); 
    end
end