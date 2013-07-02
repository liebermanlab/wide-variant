function compare_to_dnanexus(isolate_name, sorted_table_data) 
    % script for comparing to dna nexus
    % USAGE: compare_to_dnanexus('SA11-1', sorted_table_data) 

    % get dna nexus file
    filename = strcat(isolate_name, '_dnanexus.csv'); 
    variant_file = importdata(strcat('DNAnexus/', filename)); 
    
    % ---
    
    % all diverse positions called by lab pipeline
    all_positions = [sorted_table_data{:,3}]; 
    
    % ---
    
    % get which isolate this is - which site (SA11 or 12), and which isolate 
    which_site = isolate_name(3:4); 
    if strcmp(which_site, '11')
        site_num = 1;
    elseif strcmp(which_site, '12')
        site_num = 2;
    else
        fprintf('Site number not found. Try again.'); 
    end
    which_isolate = filename((strfind(filename,'-')+1):(strfind(filename,'_')-1)); 
    
    isolate_num = 24*(site_num-1) + str2num(which_isolate); 

    
    % collect lab position for each variant
    isolate_div = [sorted_table_data(:,9 + isolate_num)]; 
    isolate_div_calls = cellfun(@str2num, isolate_div); 
    isolate_lab_pos = all_positions(isolate_div_calls>0); 


    % find shared positions called between lab and dna nexus  
    shared_positions = intersect(isolate_lab_pos, variant_file); 
    num_shared_pos = numel(shared_positions); 
    num_lab_pos = numel(isolate_lab_pos); 
    lab_unique_pos = setdiff(isolate_lab_pos, variant_file); 
    
    fprintf('Isolate %s\n', isolate_name); 
    fprintf('\t%i of %i lab pipeline matched by DNAnexus\n', num_shared_pos, num_lab_pos); 
    
    % if there are lab unique positions, print them
    if numel(lab_unique_pos) ~= 0    
        fprintf('These positions are lab-pipeline unique:\n'); 
        for i = 1:length(lab_unique_pos)
            fprintf('\tPosition %i\n', lab_unique_pos(i)); 
        end
    end
end