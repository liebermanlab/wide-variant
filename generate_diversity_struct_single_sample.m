function  generate_diversity_struct_single_sample(SampleDir, SampleName, p, window_size, tempfolder)



if nargin < 5
    tempfolder = '';
end


fprintf(1,'Creating counts 3 dimensional matrix \n') ;


load([SampleDir '/diversity.mat']);


FWindows=zeros(1+window_size*2, length(p), 'single');
CWindows=zeros(1+window_size*2, length(p), 'uint16');
genome_size=size(data,2);



for j=1:length(p)
    pos=p(j);
    if pos < window_size
        st=1+window_size-pos;
        w=data(1:4,1:pos+window_size)+data(5:8,1:pos+window_size);
        [sorted, sortedpositions] = sort(w,1);
        FWindows(st+1:end,j)=(single(sorted(end,:,:))./sum(w,1))';
        FWindows(1:st,j)=0;         
        CWindows(st+1:end,j)=sum(w,1)';
        CWindows(1:st,j)=0;     
    elseif pos +window_size > genome_size
        w=data(1:4,pos-window_size:end)+data(5:8,pos-window_size:end);
        [sorted, sortedpositions] = sort(w,1);
        f=(single(sorted(end,:,:))./sum(w,1))';
        FWindows(1:length(f),j)=f;
        FWindows(length(f)+1:end,j)=0;
        CWindows(1:length(f),j)=sum(w,1)';
        CWindows(length(f)+1:end,j)=0;

    else
        %Get data
        w=data(1:4,pos-window_size:pos+window_size)+data(5:8,pos-window_size:pos+window_size);
        %Find & Store major allele frequencey
        [sorted, sortedpositions] = sort(w,1);
        FWindows(:,j)=(single(sorted(end,:,:))./sum(w,1))';
        %Store coverage
        CWindows(:,j)=sum(w,1)';

    end
end

Counts=data(:,p');

Countsi=Counts;
FWindowsi=FWindows;
CWindowsi=CWindows;

save([tempfolder '/countsatp_' SampleName], 'Countsi', 'FWindowsi', 'CWindowsi')



end

