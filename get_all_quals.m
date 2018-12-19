function all_quals_by_bp = get_all_quals(sampleinfo, genome_size) 
    % OUTPUT: where n is number of samples 
    % all_quals_by_bp: [n x GenomeLength] FQ score at each bp
    
    all_quals_by_bp = uint16(zeros(length(sampleinfo), genome_size));

    for i = 1:length(sampleinfo)    
        % get covered positions
        fprintf('\nGetting quals for isolate %i of %i\n', i, length(sampleinfo))
        
        filename = [sampleinfo(i).ExperimentFolder '/' sampleinfo(i).AlignmentFolder '/quals.mat' ];
        data = load(filename);
        all_quals_by_bp(i,:) = data.quals;
       
    end
end