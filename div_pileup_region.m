function [calls, ref1, ref2] = div_pileup_region(l, u, dir, ChrStarts)

nts=int16('ATCGatcg');

fid=fopen([char(dir) '/strain.pileup']);

line = fgetl(fid);
index=1;
calls=zeros(numel(l),4,4);
ref1=zeros(numel(l),1);
ref2=zeros(numel(l),1);

while (ischar(line) & index < numel(l)+1)
    
    l1=textscan(line,'%s','Delimiter', '\t', 'BufSize', 30000); l2=l1{1};
    chr=l2{1};
    position= double(ChrStarts(int8(chr(28))-47)+sscanf(l2{2},'%f'));
    
    if position == l(index) 

        c=int16(l2{5}); 
        %replace '.' and ',' with reference call
        ref=find(nts==int16(l2{3})); c(c==46)=nts(ref); c(c==44)=nts(ref+4);
        %replace lowercase with uppercase
        lowercase=ismember(c,nts(5:8)); c(lowercase)=c(lowercase)-32;
        
        %remove characters denoting start of reads, '^', and quality score
        startsk=find(c==94);
        c(startsk)=0;c(startsk+1)=0;
        
        %remove calls from reads ending here, '$' and actual call, which
        %precedes '$'
        endsk=find(c==36); 
        c(endsk)=0;c(endsk-1)=0;
        
        ref1(index)=ref;
        initial=c(c>0);
        
        numreads=numel(initial);
        %if all of the reads reach to the next position, we will consider
        %numreads positions. if reads end, numreads will be adjusted 
    elseif position > l(index) & position < u(index)
        
        %ignore calls from reads starting after first diverse position

        
        %remove calls from reads ending before second diverse position
        %can't just remove the first because some reads could be shorter
        %than others
        
        c=int16(l2{5}); 
        read=cumsum(ismember(c,[nts '.' ',' 'n' 'N']));
        endsr=read(c(1:find(read==numreads,1))==36);
        
        
        numreads=numreads-numel(endsr);
        initial(endsr)=[];
        

    elseif position == u(index)
        
        %process position in last line
        c=int16(l2{5}); 
        ref=find(nts==int16(l2{3})); c(c==46)=nts(ref); c(c==44)=nts(ref+4);
        lowercase=ismember(c,nts(5:8)); c(lowercase)=c(lowercase)-32;
        c(c==36)=[];
        
        final=c(1:numreads);
        
        [initialG,initialN]=ismember(initial,nts);
        [finalG,finalN]=ismember(final,nts);
        
        initialN=initialN(initialG&finalG);
        finalN=finalN(initialG&finalG);
        
        for i=1:numel(initialN)
            calls(index,initialN(i),finalN(i))=calls(index,initialN(i),finalN(i))+1;
        end
        ref2(index)=ref;    
            
        index=index+1;

        
    end
    
    line = fgetl(fid);
end



fclose(fid);

end
