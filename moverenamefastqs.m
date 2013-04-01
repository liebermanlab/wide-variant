function cmds = moverenamefastqs(IsolateTable)

global RUN_ON_CLUSTER; 
if RUN_ON_CLUSTER == 1
    mainfolder='/files/SysBio/KISHONY LAB/illumina_pipeline';
else
    mainfolder='/Volumes/sysbio/KISHONY LAB/illumina_pipeline';
end

cmds = {} ;

for i=1:length(IsolateTable)
    s = IsolateTable(i) ;
    sourcefolder = [mainfolder '/raw_data/' s.Batch '/' s.ProviderName] ; 
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