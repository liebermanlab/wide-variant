function cmds = moverenamefastqs_single_read(IsolateTable)


%Modified by TDL on 10/14/2013. Now the first column in samples.csv ('Batch') must
%contain either (1) entire path to the folder containing all fastq files or
% (2) entire path to folder containing subfolder (1 subfolder per sample)

% Modified by TDL on 10/14/2014 to accomdate the same sample being split
% across multiple batches. Requires paired end reads. Batches seperated in
% samples.csv with space

% Modified by TDL 10/2018 to enable provder name to be the name of a folder
% starting with 'ProviderName' but longer

% Modified by TDL 5/2019 to deal with unzipped starting files

cmds = {} ;

k=0;
forwardcmd=[];
reversecmd=[];

for i=1:numel(IsolateTable)
    
    iszipped=0;

    
    s = IsolateTable(i) ;
    if ~exist(s.Sample,'dir')
        mkdir(s.Sample)
    end
    
    
    if (~exist([s.Sample '/' s.Sample '.fastq'],'file')) & ~(exist([s.Sample '/sickle2050'],'dir'))
        %above line should be written more generically to test to make sure
        %alignment file isn't already there
        
        
        k=k+1;
        sourcefolders = s.Batch;
        
        
        forwardcmd=[forwardcmd ' cat '];
        %reversecmd=[reversecmd ' cat '];
        
        
        %add instructions to copy each raw file
        for j=1:numel(sourcefolders)
            
            %process folder name
            sourcefolder=sourcefolders{j};
            if sourcefolder(end)=='/'
                sourcefolder=sourcefolder(1:end-1);
            end
            
            %find fastqs
                     
            %is file in a seperate subfolder?
            f=dir([sourcefolder '/' s.ProviderName '*']); 
            if ~isempty(f) & numel(f)==1 & f.isdir > 0 %look to see if there is a folder with the provider name
                forward=dir([sourcefolder '/' f.name '/*' s.ProviderName '.fastq']);
                %reverse=dir([sourcefolder '/' f.name '/*' s.ProviderName '*2*.fastq']); 
                sourcefolder=[sourcefolder '/' f.name];
            elseif ~isempty(dir([sourcefolder '/' s.ProviderName '.fastq'])) %otherwise look in main folder
                forward=dir([sourcefolder '/' s.ProviderName '.fastq']); 
                %reverse=dir([sourcefolder '/' s.ProviderName '*2*.fastq']); 
            else %is it is .gz file?
                forward=dir([sourcefolder '/' s.ProviderName '.fastq.gz']); 
                %reverse=dir([sourcefolder '/' s.ProviderName '*2*.fastq.gz']);
                iszipped=1;
                fprintf(1,'found zipped source...')
            end
            
            
            %check
            if (numel(forward) ~= 1)
                error(['Could not find any fastq files for: ' s.ProviderName  ' in batch:' sourcefolder]);
            end
            
            forwardcmd = [forwardcmd ' ' sourcefolder '/' forward.name];
            %reversecmd = [reversecmd ' ' sourcefolder '/' reverse.name];
            
        end
        
        if k==10 %send batches of jobs to avoid sending tiny jobs to cluster
            if iszipped==0
                cmds{end+1}=[forwardcmd ' > ' s.Sample '/' s.Sample '_1.fastq'];
                %cmds{end+1}=[reversecmd ' > ' s.Sample '/' s.Sample '_2.fastq'];
            else
                cmds{end+1}=[forwardcmd ' > ' s.Sample '/' s.Sample '_1.fastq.gz; gunzip ' s.Sample '/' s.Sample '_1.fastq.gz'];
                %cmds{end+1}=[reversecmd ' > ' s.Sample '/' s.Sample '_2.fastq.gz; gunzip ' s.Sample '/' s.Sample '_2.fastq.gz'];
            end
            k=0;
            forwardcmd=[];
            %reversecmd=[];
        else
            if iszipped==0
                forwardcmd=[forwardcmd ' > ' s.Sample '/' s.Sample '_1.fastq ; '];
                %reversecmd=[reversecmd ' > ' s.Sample '/' s.Sample '_2.fastq ; '];
            else
                forwardcmd=[forwardcmd ' > ' s.Sample '/' s.Sample '_1.fastq.gz ; gunzip ' s.Sample '/' s.Sample '_1.fastq.gz; '];
                %reversecmd=[reversecmd ' > ' s.Sample '/' s.Sample '_2.fastq.gz ; gunzip ' s.Sample '/' s.Sample '_2.fastq.gz; '];
            end
      
        end
        
    end
    
    
end

if k~=0
    cmds{end+1}=forwardcmd;
    %cmds{end+1}=reversecmd;
end
