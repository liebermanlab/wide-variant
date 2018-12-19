function  generate_diversity_struct_single_sample(SampleDir, SampleName)


load('for_generate_diversity_struct_single_sample');


fprintf(1,'Creating counts 3 dimensional matrix \n') ;


load([SampleDir '/diversity.mat']);


FWindowsi=zeros(1+window_size*2, length(p), 'single');
CWindowsi=zeros(1+window_size*2, length(p), 'uint16');
genome_size=size(data,2);


Countsi=data(:,p');

if createwindows==1
    
    for j=1:length(p)
        pos=p(j);
        if pos < window_size
            st=1+window_size-pos;
            w=data(1:4,1:pos+window_size)+data(5:8,1:pos+window_size);
            [sorted, ~] = sort(w,1);
            FWindowsi(st+1:end,j)=(single(sorted(end,:,:))./sum(w,1))';
            FWindowsi(1:st,j)=0;
            CWindowsi(st+1:end,j)=sum(w,1)';
            CWindowsi(1:st,j)=0;
        elseif pos +window_size > genome_size
            w=data(1:4,pos-window_size:end)+data(5:8,pos-window_size:end);
            [sorted, ~] = sort(w,1);
            f=(single(sorted(end,:,:))./sum(w,1))';
            FWindowsi(1:length(f),j)=f;
            FWindowsi(length(f)+1:end,j)=0;
            CWindowsi(1:length(f),j)=sum(w,1)';
            CWindowsi(length(f)+1:end,j)=0;
            
        else
            %Get data
            w=data(1:4,pos-window_size:pos+window_size)+data(5:8,pos-window_size:pos+window_size);
            %Find & Store major allele frequencey
            [sorted, ~] = sort(w,1);
            FWindowsi(:,j)=(single(sorted(end,:,:))./sum(w,1))';
            %Store coverage
            CWindowsi(:,j)=sum(w,1)';
            
        end
    end
    save([TEMPORARYFOLDER '/countsatp_' SampleName], 'Countsi', 'FWindowsi', 'CWindowsi');
else
    save([TEMPORARYFOLDER '/countsatp_' SampleName], 'Countsi');
end


end

