function runtime = save_structure_parameters(logfolder, struc)
% iterate and save each field of structure to file

    runtime=datestr(now, 'yyyy-mm-dd-HH-MM-SS');
    
    % create log file with time stamp of run 
    logfilename = ['log_' runtime '.txt'];
    logfile = strcat(logfolder, '/', logfilename); 
    filehandle = fopen(logfile, 'w'); 
    
    % get all fields and values (must be numeric) 
    all_fields = fieldnames(struc); 
    for i = 1:numel(all_fields)
        fieldname = all_fields{i}; 
        class(fieldname); 
        fieldvalue = struc.(all_fields{i}); 
        class(fieldvalue);
    end
    
    fclose(filehandle); 
    
end