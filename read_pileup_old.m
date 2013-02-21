function read_pileup(fname_in,fname_out,RemoveEnds,ChrStarts,GenomeLength)

%mileup is told that the reads are in fastq format and corrects that so its
%always at Phred+33 when the mpileup comes out


fid=fopen(fname_in);

Phred_offset=33;
nts=int16('ATCGatcg');
indels=int16('+-');
numfields=26;
data=zeros(numfields,GenomeLength); %[A T C G a t c g Aq Tq Cq Gq aq tq cq gq Am Tm Cm Gm am tm cm gm D E]
% D is number of reads supporting indels
% E is number of calls at ends of a read
% Aq is the average phred qualities off all A's
% Am is the average mapping qualities of all A's


%k is the index of the character currently being read
%r is the index of the read we are currently on


mq=[]; % this holds the mapping quality of the read at position r
mq2=[]; % this holds the mapping quality of the read at position r

line = fgetl(fid);

for lklk=1:20000
%while ischar(line)

    temp=zeros(numfields,1);
    k=1; r=1; read_ends=[];
    
    lt=textscan(line,'%s','Delimiter', '\t', 'BufSize', 10000); l=lt{1};
    chr=l{1};
    position= ChrStarts((str2double(chr(28))+1))+str2double(l{2});
    ref=find(nts==int16(l{3}));
    calls=int16(l{5});
    quality=int16(l{6});
    
    

    while k < length(calls)

        
        countread=1;
        
        if calls(k)==94 %'^' %start of read, next character is a mapping quality character
            temp(end)=temp(end)+1; %records that a position was at the start or end of a read
            mq(end+1)=calls(k+1); %records mapping quality
            k=k+2; %skip mapping quality or skip mapping quality and the call
            if RemoveEnds
                countread=0;
            end
        end        
        
        if (length(calls) > k) && ~isempty(find(indels==calls(k+1),1))  %if there is an indel, we don't care if its at the front or end of a read
            k=k+2+calls(k+2)-48; %there is an additional k+1 at the bottom of this loop %the number '2' is stored as 50 in ascii
            temp(end-1)=temp(end-1)+1; %record that there was an indel
            countread=0; %don't read base that precedes an indel
        end
           
        c=calls(k);
        nt=find(nts==c);
        
        if (length(calls) > k) && calls(k+1)==36 %'$' %end of read, skip this character if RemoveEnds is on, otherwise treat as normal
            temp(end)=temp(end)+1; %records that a position was at the start or end of a read
            read_ends(end+1)=r; %will remove this position from mq at the end of this line
            k=k+1; %move one additional place in calls line so that '$' is skipped
            if RemoveEnds
                countread=0; 
            end
        end
                
        if countread==1;
            if nt
                temp(nt)=temp(nt)+1;
                temp(nt+8)=temp(nt+8)+abs(quality(r) - Phred_offset);
                temp(nt+16)=temp(nt+16)+mq(r);
            elseif c==46 %'.'
                temp(ref)=temp(ref)+1;
                temp(ref+8)=temp(ref+8)+abs(quality(r) - Phred_offset);
                temp(ref+16)=temp(ref+16)+mq(r);
            elseif c==44 %','
                temp(ref+4)=temp(ref+4)+1;
                temp(ref+12)=temp(ref+12)+abs(quality(r) - Phred_offset);
                temp(ref+20)=temp(ref+20)+mq(r);
            end
        end
        
        r=r+1; %which read (with respect to quality column and mq array)
        k=k+1; %which character
        
    end
    
    
    %calculate average scores
    for nt=1:8
        if temp(nt) > 1
            temp(nt+8)=int16(temp(nt+8)/temp(nt)); %average phred qualities
            temp(nt+16)=int16(temp(nt+16)/temp(nt)); %average mapping qualities
        end
    end
    
      
    data(:,position)=temp;
    mq(read_ends)=[];
    line = fgetl(fid);
    
end

fclose(fid);

save(fname_out,'data');

return



