function generate_diverse_positions_single_sample_snakemake( sample_path_to_variant_diversity, sample_path_to_diverse_positions, maxFQ, REF_GENOME_DIRECTORY, outgroup_boolean )

%%  Change inputs -- no vcf, just diversity.mat

%% For debugging
fprintf(1,['\n' 'Currently examining the following diversity file: ' sample_path_to_variant_diversity '\n'])
fprintf(1,['\n' 'FQ threshold: ' num2str(maxFQ) '\n'])


%% Version History
% % Tami, 2016: An earlier version of this was not readable, re-wrote
% % Arolyn, 2018.12.19: Re-written for use with snakemake pipeline, uses
% zipped files
% % Arolyn, 2018.12.20: Deletes unzipped vcf file at the end.
% % Velina modify to get diverse positions
% based on find_diverse_positions_single_sample.m and 
% generate_positions_single_sample_snakemake.m


%% For now manually add loose parameters

params = struct('minorfreqthreshold',.05, 'minreads_perstrand',1,...
    'maxreads_perstrand_percentile', 100,'minreads_perstrand_per_allele',2,...
    'min_bq',15,'min_mq',30, 'min_td', 10, 'max_td',90, 'max_sbp', 10,...
    'max_bqp', 255,'max_tdp',255, 'max_percent_ends', .90, 'max_percent_indels', .90, 'min_control_MAF', .01);

%unpack parameters
minorfreqthreshold=params.minorfreqthreshold;
%maxreads_perstrand=params.maxreads_perstrand;
minreads_perstrand=params.minreads_perstrand;
minreads_perstrand_per_allele=params.minreads_perstrand_per_allele;
min_bq=params.min_bq;
min_mq=params.min_mq;
max_td=params.max_td;
min_td=params.min_td;
max_sbp=params.max_sbp;
max_bqp=params.max_bqp;
%max_mqp=params.max_mqp;
max_tdp=params.max_tdp;
max_percent_indels=params.max_percent_indels;


num_reqs=9;

%% Get reference genome information

[ChrStarts,GenomeLength,~,ScafNames] = genomestats(REF_GENOME_DIRECTORY);

% Initialize boolean vector for positions to include as candidate SNPs that
% vary from the reference genome
p = zeros(GenomeLength,1) ;


%% For outgrouop samples only:

if outgroup_boolean
    % Then save no positions
    Positions=p2chrpos(find(p),ChrStarts);
    save([pwd '/' sample_path_to_variant_positions], 'Positions');
    return % end script, no need to look at diversity file
end


%% Get candidate SNP positions from diversity mat file

coveragethresholds=zeros(100,1);


%load data object
% -- function input must include sample_path_to_variant_diversity
load([pwd '/' sample_path_to_variant_diversity]);

%coveragethresholds -- contains cdf cutoffs .01:.01:1.0.
%only counts positions where at least 1 read aligned
%purpose of this data structure is to allow downstream processes to
%remove positions with excess coverage in a way that accounts for
%variation in coverage between samples
cov=sum(double(data(1:8,:)));
cov(cov<1)=[];


if ~isempty(cov) & sum(cov)/length(data)> 0.5

    cov=sort(cov);
    cutoffs=.01:.01:1;
    for j=1:numel(cutoffs);
        coveragethresholds(j)=cov(floor(cutoffs(j)*numel(cov)));
    end
      
    %Parse data
    NPositions=size(data,2);
    [maf, majorNT, minorNT] = div_major_allele_freq(data);
    positionsv=(1:NPositions)';
    n1=majorNT';
    n2=minorNT';
    
    minorfreq=(double(data(sub2ind(size(data),n2,positionsv)))+double(data((sub2ind(size(data),n2+4,positionsv)))))'./sum(double(data(1:8,:)));
    readsf=sum(double(data(1:4,:)));
    readsr=sum(double(data(5:8,:)));
    f1 = data(sub2ind(size(data),n1,positionsv)); %major allele counts on forward strand
    r1 = data(sub2ind(size(data),n1+4,positionsv)); %major allele counts on forward strand
    f2 = data(sub2ind(size(data),n2,positionsv)); %major allele counts on forward strand
    r2 = data(sub2ind(size(data),n2+4,positionsv)); %major allele counts on forward strand
    majorbqf=data(sub2ind(size(data),n1+8,positionsv))';
    majorbqr=data(sub2ind(size(data),n1+12,positionsv))';
    minorbqf=data(sub2ind(size(data),n2+8,positionsv))';
    minorbqr=data(sub2ind(size(data),n2+12,positionsv))';
    majormqf=data(sub2ind(size(data),n1+16,positionsv))';
    majormqr=data(sub2ind(size(data),n1+20,positionsv))';
    minormqf=data(sub2ind(size(data),n2+16,positionsv))';
    minormqr=data(sub2ind(size(data),n2+20,positionsv))';
    
    %Changed December 2012 to require forward and reverse strand qualities
    %   majorbq=(((f1.*data(sub2ind(size(data),n1+8,positionsv)))+(r1.*data(sub2ind(size(data),n1+12,positionsv))))./(f1+r1))';
    % majormq=(((f1.*data(sub2ind(size(data),n1+16,positionsv)))+(r1.*data(sub2ind(size(data),n1+20,positionsv))))./(f1+r1))';
    %  minorbq=(((f2.*data(sub2ind(size(data),n2+8,positionsv)))+(r2.*data(sub2ind(size(data),n2+12,positionsv))))./(f2+r2))';
    % minormq=(((f2.*data(sub2ind(size(data),n2+16,positionsv)))+(r2.*data(sub2ind(size(data),n2+20,positionsv))))./(f2+r2))';
    majortdF=(data(sub2ind(size(data),n1+24,positionsv)))';
    minortdF=(data(sub2ind(size(data),n2+24,positionsv)))';
    majortdR=(data(sub2ind(size(data),n1+28,positionsv)))';
    minortdR=(data(sub2ind(size(data),n2+28,positionsv)))';
    percent_indels=double(data(end,:))./double(sum(data(1:8,:))+double(data(end,:)));
    SBp=data(end-6,:);
    BQp=data(end-5,:);
    MQp=data(end-4,:);
    TDFp=data(end-3,:);
    TDRp=data(end-2,:);
    
    %Find true/false of meeting thresholds
    Tminor = minorfreq > minorfreqthreshold;
    Treads= readsf > minreads_perstrand & readsr > minreads_perstrand & (f2' > minreads_perstrand_per_allele) & (r2' > minreads_perstrand_per_allele);
    % max reads per strand removed for now in favor of thresholding later
    % & readsf < maxreads_perstrand & readsr < maxreads_perstrand ;
    
    
    %Changed December 2012 to require forward and reverse strand qualities
    Tbq= ((majorbqf > min_bq) & (minorbqf > min_bq) & (majorbqr > min_bq) & (minorbqr > min_bq));
    Tmq = ((majormqf > min_mq) & (minormqf > min_mq) & (majormqr > min_mq) & (minormqr > min_mq));
    % Tbq= (majorbq > min_bq) & (minorbq > min_bq);
    % Tmq = (majormq > min_mq) & (minormq > min_mq);
    Ttd = (majortdF > min_td) & (majortdF < max_td) & (majortdR < max_td) & (majortdR > min_td)...
        & (minortdF > min_td) & (minortdF < max_td) & (minortdR > min_td) & (minortdR < max_td);
    Tid = percent_indels < max_percent_indels;
    TSBp = SBp < max_sbp;
    TBQp = BQp < max_bqp;
    %TMQp = MQp < max_mqp;
    TTDp = (TDFp < max_tdp) & (TDRp < max_tdp);
      
    
    %Report how many positions met each requirement
%     fprintf(f,'MinorAlleleFreq: %g  \n',sum(Tminor)) ;
%     fprintf(f,'Cov: %g  \n',sum(Treads)) ;
%     fprintf(f,'minBQ: %g  \n',sum(Tbq)) ;
%     fprintf(f,'minMQ: %g  \n',sum(Tmq)) ;
%     fprintf(f,'SBp: %g  \n',sum(TSBp)) ;
%     fprintf(f,'BQp: %g  \n',sum(TBQp)) ;
%     %fprintf(f,'MQp: %g  \n',sum(TMQp)) ;
%     fprintf(f,'TDp: %g  \n',sum(TTDp)) ;
%     fprintf(f,'maxIndels: %g  \n',sum(Tid)) ;
%     fprintf(f,'acceptableTD: %g  \n',sum(Ttd)) ;
    
    
    %Records positions that met all requirements
    allreqs= Tminor + Treads + Tbq + Tmq + Ttd + Tid + TSBp + TBQp + TTDp; %TMQp
    %fprintf(f,'Max requirements met: %g  \n',num_reqs) ;
    %num_reqs=max(allreqs);
    
    fprintf(1,'Positions meeting all requirements: %g  \n',sum(allreqs==num_reqs)) ;
    
    p(allreqs==num_reqs)=1;
    
    
    % MAF control -- where does this come from?
    % seems that certain versions skip this control part 
    %Report how many positions are removed because of Isogenic control
    % skip this for now 
%     if numel(MAF_control)>1
%         good=div_single_sample_test_thresholds(data, params, MAF_control, coveragethresholds);
%         removedbycontrol=div_single_sample_test_thresholds(data, params, ones(size(MAF_control)), coveragethresholds);
%         removedbycontrol(good>0)=0;
%         fprintf(f,'Positions only removed from looseparameters because of Isogenic control: %g  \n',sum(removedbycontrol));
%     else
%         MAF_control=maf;
%     end
    
    
else    
    MAF_control=zeros(length(data),1);
end

% how many of the ~40 parameters were present in this data matrix
numfields=size(data,1);

% could figure out what to do with coveragethresholds later
coveragethresholds_sample=coveragethresholds;



% this is from other file -- remove as necessary
% % Input file
% fname_in_gz=[pwd '/' sample_path_to_variant_vcf ] % Aro_Change: print zipped filename
% gunzip(fname_in_gz); % Aro_Change: unzip
% fname_in=fname_in_gz(1:end-3) % Aro_Change: print unzipped file name
% 
% fid=fopen(fname_in,'r');
% line = fgetl(fid);
% 
% while ischar(line)
%     
%     if line(1)~='#'
%         
%         %parse line
%         lt=textscan(line,'%s','Delimiter', '\t');
%         l=lt{1};
%         
%         position_on_chr= sscanf(l{2},'%f', inf); %faster than str2double
%         position=ChrStarts(strcmp(l{1},ScafNames)) + position_on_chr;
%         
%         alt=l{5};
%         ref=l{4};
%         
%         %only consider for simple calls (not indel, not ambigious)
%         if ~isempty(alt) & ~any(alt==',') & length(alt)==length(ref) & length(ref)==1
%             
%             %find and parse quality score
%             xt=textscan(l{8},'%s','Delimiter', ';');
%             x=xt{1};
%             entrywithFQ=find(strncmp(x,'FQ',2));
%             fq=sscanf(x{entrywithFQ}((find(x{entrywithFQ}=='=')+1):end),'%f', inf);  %faster than str2double
%             
%             if int16(fq) < maxFQ; %better than maxFQ??
%                 p(position)=1;
%             end
%         end
%         
%     end
%     line = fgetl(fid);
%     
% end
% 
% fclose(fid); % close unzipped file
% delete(fname_in); % delete unzipped file


%% Save positions

Positions=p2chrpos(find(p),ChrStarts);
save([pwd '/' sample_path_to_diverse_positions], 'Positions');

end
