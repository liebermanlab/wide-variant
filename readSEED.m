function [pathways pathwaymatrix] =readSEED(filename, cstarts, glength)

%Tami Lieberman 2013


%Takes a seed annotated tsv file, ChrStarts, and Genome Length and produces
%1) A strucutre of all pathways (subsystems), containing where they were
%found, their name, and how many positions on the genome were covered
%2) A matrix of size GenomeLength X maxNumberPathwaysForGivenPositions,
%in which the elements correspond to the genomes in pathways


%Download SEED annotations from: http://pubseed.theseed.org/?page=OrganismSelect
%Organisms (select organism) -> Download - > Click here to export -> Export table

fprintf(1,'Reading in SEED annotations of pathways...');

pathways={};
pathway_names={};

pathwaymatrix=zeros(glength,1);

fid=fopen(filename);
fgetl(fid); %remove header
line=fgetl(fid);

nextpway=1;
while line~=-1
    
    delims=strfind(line,'	');
    
    ScafName=line(delims(2)+1:delims(3)-1);
    scaf=str2num(ScafName(end))+1;
    
    genestart=chrpos2index([scaf str2num(line(delims(3)+1:delims(4)-1))],cstarts);
    geneend=chrpos2index([scaf str2num(line(delims(4)+1:delims(5)-1))],cstarts);
    
    pos1=min([genestart geneend]);
    pos2=max([genestart geneend]);
    
    %disp(line)
    
    %case is empty
    
    pathwaylist=line(delims(9)+1:end);
    
    if ~strcmp(pathwaylist,'- none -	 	 ')
        pathwaydelims=[-5 strfind(pathwaylist,'; <br>') numel(pathwaylist)+1];
        for i=1:(numel(pathwaydelims)-1);
            pway=pathwaylist(pathwaydelims(i)+6:pathwaydelims(i+1)-1);
            pwaynumber=find(strcmp(pathway_names,pway));
            if isempty(pwaynumber)
                pwaynumber=nextpway;
                pathway_names{pwaynumber}=pway;
                pathways(pwaynumber).name=pway;
                pathways(pwaynumber).starts=pos1;
                pathways(pwaynumber).ends=pos2;
                nextpway=nextpway+1;
                pathways(pwaynumber).size=pos2-pos1;
            else
                pathways(pwaynumber).starts(end+1)=pos1;
                pathways(pwaynumber).ends(end+1)=pos2;
                pathways(pwaynumber).size=pathways(pwaynumber).size+pos2-pos1;
            end
            pathwaymatrix(pos1:pos2,i)=pwaynumber;
        end
        
        
    end
    line=fgetl(fid);
    %disp(size(pathwaymatrix));
end


