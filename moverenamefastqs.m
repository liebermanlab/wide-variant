function cmds = moverenamefastqs(IsolateTable)

global RUN_ON_CLUSTER; 
if RUN_ON_CLUSTER == 1
    % modified 09/03/2013 by HC
%     datafolder='/files/SysBio/KISHONY LAB/illumina_pipeline';
    datafolder = '/hms/scratch1/hattie'; 
else
    datafolder='/Volumes/sysbio/KISHONY LAB/illumina_pipeline';
end

cmds = {} ;

for i=1:length(IsolateTable)
    s = IsolateTable(i) ;
    sourcefolder = [datafolder '/raw_data/' s.Batch '/' s.ProviderName] ; 
    %fprintf(1,sourcefolder)
    if ~exist(s.Sample,'dir')
        mkdir(s.Sample)
    end
    fastqs=dir(fullfile(sourcefolder,'/*.fastq'));
    if ~exist(sourcefolder,'dir')
        error(['Could not find raw data folder called' sourcefolder])
    end    
    for j=1:numel(fastqs);
        if ~exist([s.Sample '/' s.Sample '_' num2str(j) '.fastq'],'file')
            cmds{end+1} = ['cp "' sourcefolder '/' fastqs(j).name '"  ' s.Sample '/' s.Sample '_' num2str(j) '.fastq'];
        end
    end
    
end


end