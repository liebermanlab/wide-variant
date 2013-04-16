function  [Counts, FWindows, CWindows] = generate_diversity_struct(SampleDirs, SampleNames, p, numfields, window_size, parallel, jobsubmitoptions)


global TEMPORARYFOLDER

createwindows=1;

fprintf(1,'Creating counts 3 dimensional matrix \n') ;
if (1+window_size*2 * length(p)*length(SampleDirs)) > 10^10
    fprintf(1,'Not generating FWindows and CWindows at this time -- Requires too much memory to store windows with this number of samples, window_size, and number of positions meeting loose parameters \n') ;
    createwindows=0;
    CWindows=0;
    FWindows=0;
else
    FWindows=zeros(1+window_size*2, length(p), length(SampleDirs), 'single');
    CWindows=zeros(1+window_size*2, length(p), length(SampleDirs), 'uint16');
end

Counts=zeros(numfields, length(p), length(SampleDirs),'uint16');



if parallel==1
    
    %run others
    parallel_params={};
    for i=1:numel(SampleNames)
        parallel_params{end+1}={SampleDirs{i}, SampleNames{i}, p, window_size, createwindows, TEMPORARYFOLDER};        
    end
    
    
    run_parallel_matlab_commands('generate_diversity_struct_single_sample', parallel_params, jobsubmitoptions, 1);
    
    
    %load files
    for i=1:size(SampleNames)
        load([TEMPORARYFOLDER '/countsatp_' SampleNames{i} '.mat'])
        Counts(:,:,i)=Countsi;
        if createwindows==1
            FWindows(:,:,i)=FWindowsi;
            CWindows(:,:,i)=CWindowsi;
            delete([TEMPORARYFOLDER '/countsatp_' SampleNames{i} '.mat'])
        end
    end

   
    
else

    for i=1:length(SampleDirs)


        fprintf(1,'Sample: %g  \n',i) ;
        load([SampleDirs{i} '/diversity.mat']);
        Counts(:,:,i)=data(:,p);
        genome_size=size(data,2);

        if createwindows==1
            for j=1:length(p)
                pos=p(j);
                if pos < window_size
                    st=1+window_size-pos;
                    w=data(1:4,1:pos+window_size)+data(5:8,1:pos+window_size);
                    [sorted, sortedpositions] = sort(w,1);
                    FWindows(st+1:end,j,i)=(single(sorted(end,:,:))./sum(w,1))';
                    FWindows(1:st,j,i)=0;         
                    CWindows(st+1:end,j,i)=sum(w,1)';
                    CWindows(1:st,j,i)=0;     
                elseif pos +window_size > genome_size
                    w=data(1:4,pos-window_size:end)+data(5:8,pos-window_size:end);
                    [sorted, sortedpositions] = sort(w,1);
                    f=(single(sorted(end,:,:))./sum(w,1))';
                    FWindows(1:length(f),j,i)=f;
                    FWindows(length(f)+1:end,j,i)=0;
                    CWindows(1:length(f),j,i)=sum(w,1)';
                    CWindows(length(f)+1:end,j,i)=0;
                else
                    %Get data
                    w=data(1:4,pos-window_size:pos+window_size)+data(5:8,pos-window_size:pos+window_size);
                    %Find & Store major allele frequencey
                    [sorted, sortedpositions] = sort(w,1);
                    FWindows(:,j,i)=(single(sorted(end,:,:))./sum(w,1))';
                    %Store coverage
                    CWindows(:,j,i)=sum(w,1)';

                end
            end
            
            
        end
    end
end

end

