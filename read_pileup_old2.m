function read_pileup(fname_in,fname_out,ChrStarts,GenomeLength)

%Tami Lieberman 2012

%fname_in should me the name of a pileup file output from mpileup, including
%base qualities, mapping qualities, and tail distances reported

%fname_out is what you want to save the resulting matrix as, I call it 'data.mat'

%ChrStarts is an array holding the indices in the position dimension
%corresponding to the start of a new chromsome. For B. dolosa this is [0,3413074,5589203]


% Aq is the average phred qualities of all A's
% Am is the average mapping qualities of all A's
% At is the average  tail distance of all A's
% Ps is the p value for strand bias (fishers test)
% Pb is the p value for the base qualities being the same for the two
% different types of calls (1st major, 2nd major nt, either strand) (ttest)
% Pm is the p value for the mapping qualities being the same for the two
% different types of calls (ttest)
% Pftd is the p value for the tail distantces on the forward strand
% being the same for the two different types of calls (ttest)
% Pftd is the p value for the tail distantces on the reverse strand
%being the same for the two  different types of calls (ttest)
% D is number of reads supporting indels
% E is number of calls at ends of a read


%Important constants
Phred_offset=33; %mpileup is told that the reads are in fastq format and corrects that so its always at Phred+33 when the mpileup comes out
nts=int16('ATCGatcg');
indels=int16('+-');
numfields=39;
numtests=50; %how many ttests to do in parallel

%Initialize
data=int16(zeros(numfields,GenomeLength)); %[A T C G a t c g Aq ... gq Am .... gm  At .... gt Ps Pb Pm Pftd Prtd D E]
bttests_x=nan(2000,numtests);
mttests_x=nan(2000,numtests);
fttests_x=nan(2000,numtests);
rttests_x=nan(2000,numtests);
bttests_y=nan(2000,numtests);
mttests_y=nan(2000,numtests);
fttests_y=nan(2000,numtests);
rttests_y=nan(2000,numtests);

%keep track of which test we're on for parallelizing ttests
testpositions=zeros(numtests,1);
t=0;

%holds data for each line before storing in data
temp=zeros(numfields,1);

fid=fopen(fname_in);

line = fgetl(fid);

% 2694732 was a problem line...

while ischar(line)
    
    %parse line
    lt=textscan(line,'%s','Delimiter', '\t', 'BufSize', 30000); l=lt{1};
    chr=l{1};
    %EDIT BELOW LINE
    position= ChrStarts(int8(chr(28))-47)+sscanf(l{2},'%f');
    ref=find(nts==int16(l{3}));
    calls=int16(l{5});
    bq=int16(l{6}); %base quality, BAQ corrected
    mq=int16(l{7}); %mapping quality
    td=sscanf(l{8},'%f,'); %distance from tail
    
    
    %starts of reads
    startsk=find(calls==94); %'^'
    for k=startsk
        calls(k:k+1)=-1; %remove mapping character, absolutely required because the next chracter could be $
    end
    
    %read ends
    endsk=find(calls==36); %'$'
    temp(end)=length(endsk)+length(startsk); %record how many calls were at the start/end of a read
    calls(endsk)=-1;
    
    %indels and calls from reads supporting indels
    indelk=[find(calls==43) find(calls==45)]; %'-+'
    temp(end-1)=length(indelk); %record how many indels there were
    for k=indelk
        calls(k:(k+1+(calls(k+1))-48))=-1; %don't remove base that precedes an indel
        calls(k-1)=300; %keep this number in the end, but don't coudn it towards ATCG
    end
    temp(end-1)=temp(end-1)+sum(calls==42);  %'*' indicate a missing base-- keep these in indexing, but record how you found;
    
    
    
    %replace all reference with their actual calls
    if ref
        calls(calls==46)=nts(ref); %'.'
        calls(calls==44)=nts(ref+4); %','
    end
    
    %index reads for finding scores
    simplecalls=calls((calls>0)); %simplecalls is a transformation of calls such that each calls position in simplecalls corresponds to its position in bq, mq, td
    
    
    
    %calculate how many of each type and average scores
    for nt=1:8
        if ~isempty(find(simplecalls==nts(nt),1))
            temp(nt)=sum(simplecalls==nts(nt));
            temp(nt+8)= int16(sum(bq(simplecalls==nts(nt)))/temp(nt))- Phred_offset;
            temp(nt+16)=int16(sum(mq(simplecalls==nts(nt)))/temp(nt)) - 33;
            temp(nt+24)=int16(sum(td(simplecalls==nts(nt)))/temp(nt));
        end
    end
    
    %find major and nextmajor alle
    [~, sortedpositions] = sort(temp(1:4)+temp(5:8));
    n1 = sortedpositions(end);
    n2 = sortedpositions(end-1);
    
    
    %calculate stats which require distributions during loop, strand bias
    %can de done late
    x=(simplecalls==nts(n1) | simplecalls==nts(n1+4));
    y=(simplecalls==nts(n2) | simplecalls==nts(n2+4));
    if (sum(temp(1:4)) > 20 && sum(temp(5:8))> 20 && sum(y)*200>sum(x))  %only calcualte p values of there are greater than 20 reads on each strand and MAF < .995
        t=t+1;
        %populate ttest stucture
        bttests_x(1:sum(x),t)=bq(x)-Phred_offset; bttests_y(1:sum(y),t)=bq(y)-Phred_offset;
        mttests_x(1:sum(x),t)=mq(x)-Phred_offset; mttests_y(1:sum(y),t)=mq(y)-Phred_offset;
        fttests_x(1:sum(simplecalls==nts(n1)),t)=td(simplecalls==nts(n1)); fttests_y(1:sum(simplecalls==nts(n2)),t)=td(simplecalls==nts(n2));
        rttests_x(1:sum(simplecalls==nts(n1+4)),t)=td(simplecalls==nts(n1+4)); rttests_y(1:sum(simplecalls==nts(n2+4)),t)=td(simplecalls==nts(n2+4));
        
        %record which position on chromosome these tests will
        %correspond to
        testpositions(t)=position;
        
        %calculate strand bias on the spot
        p=fexact([temp(n1) temp(n2); temp(n1+4) temp(n2+4)], 0); %strand bias
        temp(end-6)=int16(-log10(p(3)));
        
        
    end
    
    
    
    %store the data!
    data(:,position)=int16(temp);
    
    
    %if we've filled up the strucutres for tests, calculate and store in
    %data
    if t==numtests
        %calculate in parallel
        [~,bp]=ttest2(bttests_x,bttests_y);
        [~,mp]=ttest2(mttests_x,mttests_y);
        [~,fp]=ttest2(fttests_x,fttests_y);
        [~,rp]=ttest2(rttests_x,rttests_y);
        
        %record
        data(end-5,testpositions)=int16(-log10(bp));
        data(end-4,testpositions)=int16(-log10(mp));
        data(end-3,testpositions)=int16(-log10(fp));
        data(end-2,testpositions)=int16(-log10(rp));
        
        
        %clear variables
        bttests_x=nan(2000,numtests);
        mttests_x=nan(2000,numtests);
        fttests_x=nan(2000,numtests);
        rttests_x=nan(2000,numtests);
        bttests_y=nan(2000,numtests);
        mttests_y=nan(2000,numtests);
        fttests_y=nan(2000,numtests);
        rttests_y=nan(2000,numtests);
        testpositions=zeros(numtests,1);
        t=0;
    end
    
    
    %next line
    line = fgetl(fid);
    temp=zeros(numfields,1);
    
    
end

fclose(fid);

if t > 0
    %add p values for any remaining
    %calculate in parallel
    [~,bp]=ttest2(bttests_x,bttests_y);
    [~,mp]=ttest2(mttests_x,mttests_y);
    [~,fp]=ttest2(fttests_x,fttests_y);
    [~,rp]=ttest2(rttests_x,rttests_y);
    
    %record
    data(end-5,testpositions(1:t))=int16(-log10(bp(1:t)));
    data(end-4,testpositions(1:t))=int16(-log10(mp(1:t)));
    data(end-3,testpositions(1:t))=int16(-log10(fp(1:t)));
    data(end-2,testpositions(1:t))=int16(-log10(rp(1:t)));
end





save(fname_out,'data');

return

