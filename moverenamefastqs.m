function cmds = moverenamefastqs(IsolateTable)

cmds = {} ;

for i=1:length(IsolateTable)
    s = IsolateTable(i) ;
    sourcefolder = ['../raw_data/' s.Batch '/' s.ProviderName '/'] ; 
    %fprintf(1,sourcefolder)
    if ~exist(s.Sample,'dir')
        mkdir(s.Sample)
    end
    fastqs=dir(fullfile(sourcefolder,'*.fastq'));
    for j=1:numel(fastqs);
        if ~exist([s.Sample '/' s.Sample '_' num2str(j) '.fastq'],'file')
            cmds{end+1} = ['cp ' sourcefolder fastqs(j).name '  ' s.Sample '/' s.Sample '_' num2str(j) '.fastq'];
        end
    end
    
end


end