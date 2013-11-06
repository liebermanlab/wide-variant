function analysis_ambiguous_N_tree(Nthreshold, Calls, p, SampleNames, MutQual, qual_0)
    % params
    mainfolder = '/Volumes/sysbio/KISHONY LAB/illumina_pipeline';
    RefGenomeNameUnedited = 'Smaltophilia_K279';
    num_trees = 10;  
    
    num_pos = 1:size(Calls,1); 
    quality_positions = num_pos(MutQual>qual_0);
    quality_positions_chromosomal = p(quality_positions); 
    Callsgood = Calls(quality_positions,:);
    isN = get_ambiguous_calls(Callsgood); 
    
    % Creates tree input file based on N threshold  
    lowNpositions = sum(isN,2)<=Nthreshold; 
    calls_for_tree = Callsgood(lowNpositions,:);
    
    % ADD REFERENCE AT THESE POSITIONS
    [outgroup_nts, ~] = extract_outgroup_mutation_positions(mainfolder, RefGenomeNameUnedited, quality_positions_chromosomal); 
    outgroup_nts = outgroup_nts(lowNpositions); 
    calls_for_tree = [outgroup_nts', calls_for_tree]; 
    outgroup_num = size(Calls,2)+1; 
    
    TreeSampleNames{1} = 'Outgroup'; 
    TreeSampleNames(2:length(SampleNames)+1) = SampleNames; 
    
    %write input file
    filenamebit = ['Nthreshold_' num2str(Nthreshold)]; 
    
    generate_phylip_input(calls_for_tree, TreeSampleNames, [filenamebit '_infile.txt'])

    %write option file
    fid = fopen([filenamebit '_optionfile.txt'],'w');
    fprintf(fid, [filenamebit '_infile.txt\n']);
    fprintf(fid, 'f\n');
    fprintf(fid, [filenamebit '_out.txt\n']);
    fprintf(fid, 'v\n'); 
    fprintf(fid, '%i\n', num_trees); 
    fprintf(fid, 'o\n'); 
    fprintf(fid, '%i\n', outgroup_num); 
    fprintf(fid, 'y\n');
    fprintf(fid, 'f\n'); 
    fprintf(fid, [filenamebit '_' num2str(num_trees) 'trees.tree\n']); 
    fclose(fid);

    %run
    fid = fopen('temp.sh','w');
    fprintf(fid, ['! ./dnapars <' filenamebit '_optionfile.txt >' filenamebit '_outfile.txt\n']);
    fclose(fid);
    
    fprintf('\nGenerating %i trees for N threshold %i\n', num_trees, Nthreshold); 
    ! chmod +x temp.sh
    ! ./temp.sh
    fprintf('\nFinished.\n'); 

end