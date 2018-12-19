function vcf_to_quals_16S(sample_directory,REF_GENOME_DIRECTORY)


%% Version history

% Tami Lieberman 2012, modified 2016

% Given a vcf file with one file per line, grabs the FQ score for each
% position. Ignores lines corresponding to indels

% 2018.10.16, Arolyn: New function to grab strain_16S.vcf rather than
% strain.vcf, in order to get 16S gene info; write quals_16S.mat rather
% than quals.mat


%%

fname_in=[sample_directory '/strain_16S.vcf'];

if nargin < 2
    load for_matlab_scripts
    path(SCRIPTSDIRECTORY,path);
end

[ChrStarts,GenomeLength,~,ScafNames]=genomestats(REF_GENOME_DIRECTORY);


quals=zeros(1,GenomeLength,'int16');


i=1;
position=0;

fid=fopen(fname_in,'r');
line = fgetl(fid);


while ischar(line)
    
    if line(1)~='#' 
        
        if ~mod(i,50000), fprintf(1,'.'); end
        
        %parse line
        lt=textscan(line,'%s','Delimiter', '\t');
        l=lt{1};
        
        chr= find(strcmp(l{1},ScafNames));
        position_on_chr= sscanf(l{2},'%f', inf); %faster than str2double
        position=ChrStarts(chr) + position_on_chr;
                
        alt=l{5};
        ref=l{4};
        
        %only store quality for simple calls (not indel, not ambigious)
        if ~isempty(alt) & ~any(alt==',') & length(alt)==length(ref) & length(ref)==1
            
            %find and parse quality score
            xt=textscan(l{8},'%s','Delimiter', ';');
            x=xt{1};
            
            
            entrywithFQ=find(strncmp(x,'FQ',2));
            fq=sscanf(x{entrywithFQ}((find(x{entrywithFQ}=='=')+1):end),'%f', inf);  %faster than str2double
            
            if fq < quals(position); %if already a position with a stronger FQ here, don't include this. (More negative is stronger)
                quals(position)= int16(fq);
            end
        end
        
        i=i+1;
    end
    line = fgetl(fid);
    
end

fclose(fid);
save([sample_directory '/quals_16S.mat'], 'quals');

end

