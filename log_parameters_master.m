% add promoter size, qual_0 to log

log_parameters = STRICT_PARAMETERS; 
log_parameters.('promotersize') = promotersize; 
log_parameters.('qual_o') = qual_0; 

logfolder = 'log_analysis/'; 
if exist(logfolder, 'dir') ~= 7
    mkdir(logfolder)
end

save_structure_parameters(logfolder, log_parameters); 
