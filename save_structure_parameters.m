function save_structure_parameters(logfolder, struc)
% iterate and save each field of structure to file

    % create log file with time stamp of run 
    logfilename = ['log_' datestr(now, 'yyyy-mm-dd-HH-MM-SS') '.txt'];
    logfile = strcat(logfolder, '/', logfilename); 
    filehandle = fopen(logfile, 'w'); 
    
    % get all fields and values (must be numeric) 
    all_fields = fieldnames(struc); 
    for i = 1:numel(all_fields)
        fieldname = all_fields{i}; 
        class(fieldname); 
        fieldvalue = struc.(all_fields{i}); 
        class(fieldvalue)
        fprintf(filehandle, '%s,%d\n', fieldname, fieldvalue)
    end
    
    fclose(filehandle); 
    
end