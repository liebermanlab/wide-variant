function pileup_to_diversity_matrix_light(dir)

if nargin > 0
    cd(dir)
end

fprintf(1,pwd)

%Tami Lieberman 2012

%edited Seungsoo Kim January 2013 to account for indels of size >= 10
%data type changed to uint8

%edited by Tami Lieberman 2014 to run from commandline


%% Important constants
Phred_offset=33; %mpileup is told that the reads are in fastq format and corrects that so its always at Phred+33 when the mpileup comes out
nts=int16('ATCGatcg');
indels=int16('+-');
numfields=8;

%% Find out which reference genome we are in

load alignment_info
refgenome=ae.Genome;

fname_in='strain.pileup';
fname_out='diversity.mat';

%% Initialize
[ChrStarts,GenomeLength,ChromosomeIndicator,ScafNames]=genomestats(['/scratch/mit_lieberman/Reference_Genomes/' refgenome]); %RUN_ON_CLUSTER
data=zeros(numfields,GenomeLength,'uint16'); %[A T C G a t c g Aq ... gq Am .... gm  At .... gt Ps Pb Pm Pftd Prtd E D]

fid=fopen(fname_in);

line = fgetl(fid);
linec = 0; 

while ischar(line)
	%fprintf(1, [line '\n']);     
    %parse line
    lt=textscan(line,'%s','Delimiter', '\t'); l=lt{1};
    chr=l{1}; 
    
    if numel(ChrStarts)==1
        position= sscanf(l{2},'%f', inf);
    else
        if sum(ismember(ScafNames,chr))==0
            error('Scaffold name in pileup file not found in reference')
        end
        %e.g. define position by looking at 28th character in string defining B. dolosa chrosome
        position= double(ChrStarts(ismember(ScafNames,chr))+sscanf(l{2},'%f'));
    end
  
    
    ref=find(nts==int16(l{3}));
    if ref > 4
        ref = ref - 4;
    end

    calls=int16(l{5});
    
    %starts of reads
    startsk=find(calls==94); %'^'
    for k=startsk
        calls(k:k+1)=-1; %remove mapping character, absolutely required because the next chracter could be $
    end
    
    %read ends
    endsk=find(calls==36); %'$'
    calls(endsk)=-1;
    
    %indels and calls from reads supporting indels
    indelk=[find(calls==43) find(calls==45)]; %'-+'
    for k=indelk
        %SK: if indel is of size >= 10, take next two characters
        if calls(k+2) >= 48 && calls(k+2) < 58
            indelsize=str2double(char(calls(k+1:k+2)));
            indeld=2;
        else
            indelsize=double(calls(k+1)-48);
            indeld=1;
        end
        %remove indel calls from counting
        calls(k:(k+indeld+indelsize))=-1; %don't remove base that precedes an indel
        calls(k-1)=300; %keep this number in the end, but don't coudn it towards ATCG
    end
    %ignore '*', as these deletion markers won't be counted towards score
    %qualities, etc and the presence of a deletion is recorded by the
    %upstream '-' indicator
    
    
    
    %replace all reference with their actual calls
    if ref
        calls(calls==46)=nts(ref); %'.'
        %(1,[num2str(ref) ',' num2str(position) '\n']); %for error tracking
        calls(calls==44)=nts(ref+4); %','
    end
    
    %index reads for finding scores
    simplecalls=calls((calls>0)); %simplecalls is a transformation of calls such that each calls position in simplecalls corresponds to its position in bq, mq, td
    
    
    
    %calculate how many of each type and average scores
    for nt=1:8
        if ~isempty(find(simplecalls==nts(nt),1))
            data(nt,position)=sum(simplecalls==nts(nt));
        end
    end
    
    line = fgetl(fid);
	linec = linec+1;     
    
end

fclose(fid);


coverage=sum(data(1:8,:));

save('coverage','coverage');

save(fname_out,'data');

return

