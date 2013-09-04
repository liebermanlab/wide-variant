function [CStarts, GLength, CIndicator]= genomestats(genome, RUN_ON_CLUSTER)

if RUN_ON_CLUSTER==1
   fr = fastaread(['/home/hc168/Reference_Genomes/' genome '/genome.fasta']) ;
else  
    fr = fastaread(['/Volumes/sysbio/KISHONY LAB/illumina_pipeline/Reference_Genomes/' genome '/genome.fasta']) ;
end

GLength=0;
CStarts=[];

ScafNames = {fr.Header} ;
for i=1:length(ScafNames)
    f=find(ScafNames{i}==' ',1) ;
    ScafNames{i} = ScafNames{i}(1:f-1) ;
    CStarts(end+1)=GLength;
    GLength=GLength+numel(fr(i).Sequence);
end

CIndicator = find(fr(1).Header=='.',1) - 1; %Assumes fewer than 10 chromosomes
