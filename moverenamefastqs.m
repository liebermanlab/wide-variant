function cmds = moverenamefastqs(IsolateTable)

%Modified by TDL on 10/14/2013. Now the first column in samples.csv ('Batch') must
%contain either (1) entire path to the folder containing all fastq files or
% (2) entire path to folder containing subfolder (1 subfolder per sample)

cmds = {} ;

for i=1:length(IsolateTable)
    
    s = IsolateTable(i) ;
    
    sourcefolder = s.Batch;
    if sourcefolder(end)=='/'
        sourcefolder=sourcefolder(1:end-1);
    end
    
    if ~exist(s.Sample,'dir')
        mkdir(s.Sample)
    end
    
    %find fastqs
    if exist([sourcefolder '/' s.ProviderName],'dir')
        fastqs=dir([sourcefolder '/' s.ProviderName '/*.fastq']); % HC fixed bug
        sourcefolder=[sourcefolder '/' s.ProviderName];
    else
        fastqs=dir([sourcefolder '/' s.ProviderName '/*.fastq']); % HC fixed bug
    end
    %copys
    if numel(fastqs) < 1
        error(['Could not find fastq file(s) for: ' s.ProviderName]);
    end
    
    for j=1:numel(fastqs);
        if ~exist([s.Sample '/' s.Sample '_' num2str(j) '.fastq'],'file')
            cmds{end+1} = ['cp "' sourcefolder '/' fastqs(j).name '"  ' s.Sample '/' s.Sample '_' num2str(j) '.fastq'];
        end
    end
    
end


end
