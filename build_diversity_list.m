%% Important variables

%Control is first sample in case


path('../scripts',path)
deep_consensus_directory='../../../Illumina_pipeline/Case_2012_11_deep_seq_from_sputum_calls2';
isolate_directory='../../../Illumina_pipeline/Case_2012_11_isolates_from_sputum2';
pt2NatGenfile='../../../Illumina_pipeline/Case_2012_Natgen_Pt_J/patientJmutations.mat';
filen='diversity.mat';
postfix='12_12_20';

comparePrevStudy=1;
usefixedpos=1; fixedscore=282;
useisolates=1;
processed=0;
showGenomicDistributions=0;
showScattersforIsogenicControl=0;

ChrStarts=[0,3413073,5589202];
GenomeLength=6420400;
window_size=500;
promotersize=150;
NTs='ATCG';

scrsz = get(0,'ScreenSize');

samples=[6 2 5 7 4];
nonmutatorsamples=[6 2 5 7];
mutator=4;
poscontrolsample=6;
sampleclosesttoroot=5;
comparisonpt=2; %comparison to previous isolates;
noncomparisionsamples=[6 5 7 4];


referenceisancestor=1;


%figureBlue = [.14 .43 .77];
%tamisgreen=[0 0 .77];

tamismagenta=[.85 .31 .85];
tamisgreen=[145 220 83]/256;


prevmutated=[1161 1160 4180 2650 4154 93 2781 2180 3694 4506 220 2317 2321 1250 856 3003 00179];




%% Set filtering parameters

%read_pileup cuts off at .995 frequency
%All parameters are strict > and < 


%loose parameters doesn't have an upper coverage thershold yet
loose_parameters=struct('minorfreqthreshold',.02, 'minreads_perstrand',30,...
    'minreads_perstrand_per_allele',2, 'min_bq',15,'min_mq' ,30, 'min_td', 0,...
    'max_td',50, 'max_sbp', 10,'max_bqp', 255,'max_tdp',255, 'max_percent_indels', .20, 'min_control_MAF', .01);


%maxreads_perstrand_percentile is which threshold in list .01:.01:1 ...
%    e.g. 98 is 98 percentile of covered positions

strict_parameters=struct('minorfreqthreshold',.03, 'maxreads_perstrand_percentile',99,...
    'minreads_perstrand',30, 'minreads_perstrand_per_allele',2,'min_bq',19,'min_mq', 33, 'min_td', 7,...
    'max_td',42, 'max_sbp', 3,'max_percent_indels', .20, 'min_control_MAF', .98, ...
    'max_bqp', 200,'max_tdp',200, 'max_percent_ends', 1);

% last 3 are essentially not used
%  
% test_parameters=struct('minorfreqthreshold',.03, 'maxreads_perstrand_percentile',99,...
%     'minreads_perstrand',30, 'minreads_perstrand_per_allele',2,'min_bq',19,'min_mq', 33, 'min_td', 7,...
%     'max_td',42, 'max_sbp', 3,'max_percent_indels', .20, 'min_control_MAF', .01, ...
%     'max_bqp', 200,'max_tdp',200, 'max_percent_ends', 1);
% 
% test2_parameters=struct('minorfreqthreshold',.03, 'maxreads_perstrand_percentile',99,...
%     'minreads_perstrand',30, 'minreads_perstrand_per_allele',2,'min_bq',10,'min_mq', 33, 'min_td', 7,...
%     'max_td',42, 'max_sbp', 3,'max_percent_indels', .20, 'min_control_MAF', .01, ...
%     'max_bqp', 200,'max_tdp',200, 'max_percent_ends', 1);
% goodpos=div_test_thresholds(counts,strict_parameters, coveragethresholds);
% testpos=div_test_thresholds(counts,test_parameters,coveragethresholds);
% goodpos(fixedpos>0)=0; 
% testpos(fixedpos>0)=0; 
% sum(testpos-goodpos)
% x=find(testpos(:,6)-goodpos(:,6)~=0)
% %intersect_iMAF(testpos(:,6)-goodpos(:,6)~=0)
% p2chrpos(p(sum(testpos,4)-sum(goodpos,4)~=0),ChrStarts)
% testpos2=div_test_thresholds(counts,test2_parameters,coveragethresholds);
% testpos2(fixedpos>0)=0; 


 
%% Build position list, display histrograms of each quality score

switch 0
    case 0

        % Get SampleNames and ScafNames
        SampleNames = read_sample_names ;
        NSample = length(SampleNames) ;
        RefGenome = {} ;
        SampleDirs = {} ;
        
        for i=1:NSample
            SampleDirs{i} = ['../' SampleNames(i).ExperimentFolder '/' SampleNames(i).Sample '/' SampleNames(i).AlignmentFolder ] ;
            ainfo = load([SampleDirs{i} '/alignment_info']) ;
            RefGenome{i} = ainfo.ae.Genome ;
        end

        RefGenome = unique(RefGenome) ;
        if length(RefGenome)>1
            error('Must compare samples aligned to the same reference genome')
        end
        
        RefGenome = RefGenome{1} ;
       
        
        fr = fastaread(['../Reference_Genomes/' RefGenome '/genome.fasta']) ;
        
        ScafNames = {fr.Header} ;
        for i=1:length(ScafNames)
            f=find(ScafNames{i}==' ',1) ;
            ScafNames{i} = ScafNames{i}(1:f-1) ;
        end
        
        
        %Find diverse positions
        [dp, positions, numfields, coveragethresholds] = find_diverse_positions(filen, loose_parameters, strict_parameters, SampleDirs, SampleNames, ChrStarts, showGenomicDistributions, showScattersforIsogenicControl);
       
        %Save things
        savename=['Positions_' postfix];
        if useisolates==1
          % Load in isolate information to save alongsize
          load([isolate_directory '/isolatesfinal'])
          save(savename, 'dp', 'positions', 'numfields', 'coveragethresholds', 'SampleNames', 'isolateCalls', 'isolateNames', 'isolatesCoverage', 'isolatePositions', 'ScafNames', 'ainfo', 'SampleDirs', 'NSample', 'RefGenome')
        else
           save(savename, 'dp', 'positions', 'numfields', 'coveragethresholds', 'SampleNames', 'ScafNames', 'ainfo', 'SampleDirs', 'NSample', 'RefGenome')
        end
        
        
    case 1
        savename=['Positions_' postfix];
        load(savename)
end


fprintf(1,'Found %g diverse positions that meet specified criteria \n',length(dp)) ;


 



%% Load isolates

if useisolates==1
    
    [isolateMAF, isolateMinorAF,isolateMajorNT, ~]=strains_major_allele_freq(isolateCalls);

    %sort isolate positions and related structures
    ip=chrpos2index(isolatePositions', ChrStarts);
    [ip, sortedpositions] = sort(ip);
    iMAF=isolateMAF(sortedpositions);
    iMinorAF= isolateMinorAF(sortedpositions);
    iMajorNT=isolateMajorNT(sortedpositions);
    isolatePositions=isolatePositions(sortedpositions,:);
    isolateCalls=isolateCalls(sortedpositions,:);

    
else
    ip=[];
end





%% Load deep sequencing consensus
%fixed means that the allele went to fixation in that population

if usefixedpos==1
    load([deep_consensus_directory '/consensusfinal'])
    [cp, sortedpositions] =sort(chrpos2index(consensusPositions', ChrStarts));
    consensusScore=consensusScore(sortedpositions,:);
else
    cp=[];
end

%% Load info about Patient J/ P2 from Nature Genetics study

if comparePrevStudy==1
    load(pt2NatGenfile)
    NGfixedp=chrpos2index(fixedPositionsNatgen', ChrStarts);
    NGpolymorphicp=chrpos2index(polymorphicPositionsNatgen',ChrStarts);
else
    NGfixedp=[];NGpolymorphicp=[];
end



%% Combine p from isolates and deep and deep consensus


dp=unique([dp; cp]);
allp=unique([dp; ip; NGpolymorphicp; NGfixedp;]);
p=sort(allp);
positions=p2chrpos(p,ChrStarts);
consensusCalled=zeros(numel(p),size(consensusCalls,2));
consensusCalled(ismember(p,cp),:)=consensusScore>=fixedscore;
otherp=unique([NGpolymorphicp; NGfixedp]);

%% Build counts list, annotate mutations

tic;

switch processed
    case 0
        [counts, fwindows, cwindows] = generate_diversity_struct(filen, SampleDirs, p, numfields, window_size) ;
        [geneloc,  cds, mutations, sequences] = annotate_mutations_auto_gb(positions,ScafNames,RefGenome) ;
        
        savename1=['Counts_annotation_' postfix];
        savename2=['Windows' postfix];
        save(savename1,'counts', 'mutations', 'geneloc', 'cds', 'sequences')
        save(savename2, 'fwindows', 'cwindows')
    case 1
        savename1=['Counts_annotation_' postfix];
        savename2=['Windows' postfix];
        load(savename1)
        load(savename2)
end
toc;





%% Make isolate windows
if useisolates==1
    [intersect_iMAF, intersect_iMinorAF, intersect_iMajorNT, isolate_windows]=isolate_MAF_in_region(p, ip,iMAF, iMinorAF, iMajorNT, ChrStarts, window_size);


    %This inverts isolate_windows based on maf of deep sequencing from S2,
    %isolate windows is only used for graphing (changes to -1 for plotting)
    %We do not change the the record of the isolate major allele or
    %anything other than isolate_windows
    [~, maNT, ~] = div_major_allele_freq(counts); %also done later
    goodpos=div_test_thresholds(counts,strict_parameters, coveragethresholds); %also done later
    reverseMajors=find(maNT(:,poscontrolsample)~=intersect_iMajorNT & goodpos(:,poscontrolsample)>0 & intersect_iMajorNT~=0);
    for i=1:numel(reverseMajors)
        tempm=cell2mat(isolate_windows(reverseMajors(i)));
        tempm(2,tempm(1,:)==p(reverseMajors(i)))=-1*tempm(2,tempm(1,:)==p(reverseMajors(i))); %reverse position that matches the center
        isolate_windows(reverseMajors(i))={tempm};
    end
    
end




%% Table

if useisolates==1
    [combined_annotation_goodpos, combined_annotation_all, updated_annotation_diverse,  updated_annotation_isolates] = div_clickable_table_isolates(mutations, dp, ip, otherp, isolateCalls, counts,  fwindows, cwindows, isolate_windows, strict_parameters, coveragethresholds, RefGenome, ScafNames, SampleNames, ChrStarts, promotersize);
else
    [combined_annotation_goodpos, combined_annotation_all] = div_clickable_table(mutations, p, counts,  fwindows, strict_parameters, RefGenome, ScafNames, SampleNames, ChrStarts, promotersize);
end

 

genes={combined_annotation_all.gene};
types=[combined_annotation_all.type];
itypes=[updated_annotation_isolates.type];
typesmatrix=repmat(types,7,1)';
%igenes={updated_annotation_isolates.gene),1)';

genes2=genes; genes2(cellfun(@numel,genes)<2)={'BDAG_00000'};
gm=cell2mat(genes2');
gn=str2num(gm(:,6:end));
gnmat2=repmat(gn,1,7);


%[combined_annotation_S, ~, ~, ~] =div_clickable_table_only_S(mutations, dp, ip, isolateCalls, counts,  fwindows, isolate_windows, strict_parameters, RefGenome, ScafNames, SampleNames, ChrStarts);






%% Calculate important things, directionality, goodpos, fixedpos

[maf, maNT, minorNT] = div_major_allele_freq(counts);
minorAF=div_minor_allele_freq(counts);
[~,refnt]=ismember([combined_annotation_all.ref],NTs);

goodpos=div_test_thresholds(counts,strict_parameters, coveragethresholds);


if referenceisancestor
    ancnt=refnt';
    ancnt(p==2849264)=1; %glycosyltransferase
    % if reference isn't found in any of our samples, use mode
    for i=1:numel(allp)
        if isempty(find(maNT(i,samples)==refnt(i),1))
            ancnt(i)=mode(maNT(i,samples));
            %if its a 3/2 split, go with what pt 3 has (based on phylogeny), otherwise use mode
            
            if sum(maNT(i,samples)==ancnt(i))==3
                ancnt(i)=maNT(i,sampleclosesttoroot);
            end
        end
    end
else
    ancnt=mode([maNT(:,samples) refnt'],2);
end


% find positions samples in which alternative allele is fixed
% force fixedpos called if it was called in another strain (cpm >0) and mutation allele frequency is greater than .95

refntm=repmat(refnt',1,7);
ancntm=repmat(ancnt,1,7);
cpi=(ismember(p,cp)); cpm=repmat(cpi,1,7);
fixedgenes=cell2mat(genes(cpi' & cellfun(@numel,genes) > 0)');
fixedgenesN=str2num(fixedgenes(:,6:end));
fixedpos=cpm & ancntm~=maNT & maf > .95; 
%any 'goodpos' that are also 'fixedpos' are now no longer diverse
goodpos(fixedpos>0)=0; 


%find allele frequency of mutant (non-reference!!!) allele
mutAF=minorAF;
mutAF(ancntm~=maNT)=maf(ancntm~=maNT);
intersect_iMutAF= intersect_iMinorAF;
intersect_iMutAF((ancnt~=intersect_iMajorNT) & intersect_iMajorNT >0)=  intersect_iMAF((refnt'~=intersect_iMajorNT) & intersect_iMajorNT >0);
intersect_iMutAF(fixedpos(:,poscontrolsample)>0)=1; %this step required because the isolate pipeline doesn't currently keep fixed positions in isolates in this pipe --  removes nonpolymorphic positions before saving



%% Shared vs private fixed mutations


uniquefixedpos=zeros(1,numel(SampleDirs));
sharedfixedpos=zeros(1,numel(SampleDirs));
privatefixedmutation=zeros(size(fixedpos));
for i=samples
    othersamples=samples; othersamples(othersamples==i)=[];
    fpos=find(fixedpos(:,i)>0);
    for k=1:numel(fpos)
        j=fpos(k);
        if sum(fixedpos(j,othersamples))>0 & ~isempty(maNT(j,fixedpos(j,othersamples)>0)==maNT(j,i))
            sharedfixedpos(i)=sharedfixedpos(i)+1;
        else
            uniquefixedpos(i)=uniquefixedpos(i)+1;
            privatefixedmutation(j,i)=1;
        end
    end
end



figure(472); clf; hold on; 
h=bar([uniquefixedpos(nonmutatorsamples) 28; sharedfixedpos(samples)]','grouped') %supress mutator sample because it has 
%set(gca,'Xticklabel',{SampleNames(samples).Sample})
colormap([rgb('Grey'); rgb('Black')])
%set(h,'EdgeColor', 'k')
legend({'Patient-specific fixed mutations', 'Fixed mutations shared with other patients'})
grid(gca)
set(gca,'Xtick',1:5)
set(gca,'Xticklabel',{'P1', 'P2', 'P3', 'P4', 'P5'})
ylabel('Fixed mutations')

set(gca,'Ytick',[0:5:20 28])
set(gca,'Yticklabel',[0:5:20 uniquefixedpos(mutator)])

set(gca, 'GridLineStyle', ':');


%phylip header
disp(['   ' num2str(numel(samples)+1) '   ' num2str(sum(sum(fixedpos(:,samples),2)>0))])
k=1;
for i=samples
    fcalls=NTs(ancnt);
    fcalls(fixedpos(:,i)>0)=NTs(maNT(fixedpos(:,i)>0,i));
    fcalls=fcalls(sum(fixedpos(:,samples),2)>0);
        
    disp(['P' num2str(k) '        ' fcalls]);
    k=k+1;
end
disp(['ref       ' NTs(refnt(sum(fixedpos(:,samples),2)>0))])







%% Compare NatGen and new results for patient 2
%positions are fixed if they are fixed relative to P1 (S2)


%this chunk of code was copied into div_natgen_diversity_comparison for
%storage and notes 


%also see Illumina_pipeline/Case_2012_Natgen_Pt_J/ana_master.m

%% Specificity analysis for allele frequency using multiple colonies

%A chunk of code was put into a function
%div_exhaustive_sensitivity_specificity
%Frequency ROCAUC

minorfreqthreshold=0:.0025:1; 
test_parameters=strict_parameters;
x=[]; y=[];
for i=minorfreqthreshold
    test_parameters.minorfreqthreshold=i;
    
    [spec, sens, normalized_spec, expected_spec]=div_specificity_sensitivity(intersect_iMutAF, mutAF(:,poscontrolsample), maf(:,1), counts(:,:,poscontrolsample), test_parameters, coveragethresholds(:,poscontrolsample), maNT(:,poscontrolsample)', minorNT(:,poscontrolsample)');


    x(end+1)=1-spec;
    y(end+1)=sens;
end
 


%change number of colonies required

numcolonies=1:1:29;
x3=[]; y3=[];
for i=numcolonies
    
    %repeat only with polymorphisms found in at least 2 colonies
    intersect_iMutAF2=intersect_iMutAF; intersect_iMutAF2(intersect_iMutAF2<i/29)=0;
    [spec, sens, normalized_spec, expected_spec]=div_specificity_sensitivity(intersect_iMutAF2, mutAF(:,poscontrolsample), maf(:,1), counts(:,:,poscontrolsample), strict_parameters, coveragethresholds(:,poscontrolsample), maNT(:,poscontrolsample)', minorNT(:,poscontrolsample)');
    % plot(1-normalized_specificity,sensitivity,'r.','MarkerSize', 5, 'ButtonDownFcn',{@dispparamters, test_parameters});
   % plot(1-specificity,sensitivity,'k.','MarkerSize', 5, 'ButtonDownFcn',{@dispparamters, test_parameters});
    x3(end+1)=1-spec;
    y3(end+1)=sens;
end


figure(9); clf; hold on;
plot(-0.1:.1:1,-0.1:.1:1,':k')
plot(x,y, 'k', 'LineWidth', 2, 'Color', rgb('Gray'))
plot(x3,y3, 'k--', 'LineWidth', 2, 'Color', rgb('Sienna'))
for i=[.15 .10 .05 .03 .02 .01]
    j=find(minorfreqthreshold==i);
    plot(x(j),y(j),'ko','MarkerFaceColor', rgb('Gray'),'MarkerSize', 6)
end   
for i=[1 2 3 5 ]
    j=find(numcolonies==i);
    plot(x3(j),y3(j),'ko','MarkerFaceColor', rgb('Brown'),'MarkerSize', 6)
end   

plot(x3(1),y3(1),'kd', 'Color', 'k', 'MarkerFaceColor','k', 'MarkerSize', 8)
axis([-.01 1 -.01 1])


%% Some Summarys
figure;
bar(1:numel(samples),[sum(fixedpos(:,samples)); sum(goodpos(:,samples))]')
legend({'Fixed', 'Polymorphic'})
ylabel('Number of SNPs')
set(gca,'Xtick',1:numel(samples))
set(gca,'Xticklabel',{SampleNames(samples).Sample})


figure(3007); clf; hold on;
bar(1:numel(samples),[sum(goodpos(:,samples))]',.5)
colormap([rgb('Black')])
%legend({'Fixed', 'Polymorphic'})
ylabel('Number of SNPs')
set(gca,'Xtick',1:numel(samples))
set(gca,'Xticklabel',{'P1', 'P2', 'P3', 'P4', 'P5'})
%set(gca,'Xticklabel',{SampleNames(samples).Sample})
grid(gca)
set(gca, 'GridLineStyle', ':');


%% Make a convenient data structure for accessing size of genes

%get coding sequence numbers
allcds=[cds{:}];
allgenenames=char({allcds.gene}');
allgenenums=str2num(allgenenames(:,6:end));

%get sizes
%first order cds
[allgenenums,codingsequencepostions]=sort(allgenenums);
allcds=allcds(codingsequencepostions);

for i=find(cellfun(@isempty,{allcds(:).loc2})>0)
    allcds(i).loc2=allcds(i).loc1+1000; %make arbitrary sizes for others
end

genesizes=zeros(1,max(allgenenums));
genesizes(allgenenums)=[allcds.loc2]-[allcds.loc1];



%% Display clickable scatter for data analysis of MAF versus MAF of control

controlFreq=maf(:,1);
controlNT=maNT(:,1);


for i=samples
    div_clickable_scatter_sigcolor(maf(:,1), maf(:,i), 'Control MAF', ['MAF in ' SampleNames(i).Sample], i, strict_parameters, coveragethresholds, counts, fwindows, isolate_windows, cwindows, positions, mutations, RefGenome, ScafNames,  ChrStarts, SampleNames);
end



%% Scatter plot of comparision to individual isolates, for sputum #2 (sample 6)





bothmethods=(intersect_iMinorAF > 0 | goodpos(:,poscontrolsample)>0| fixedpos(:,poscontrolsample)>0);

bothmethods03=(intersect_iMinorAF > 0.03 | goodpos(:,poscontrolsample)>0 | fixedpos(:,poscontrolsample)>0);


div_clickable_scatter(intersect_iMutAF, mutAF(:,poscontrolsample),  'Sputum#2 isolates mutation frequency','Sputum#2 deep mutation frequency',6, bothmethods, strict_parameters, coveragethresholds, counts, fwindows, isolate_windows, cwindows, positions, mutations, RefGenome, ScafNames, ChrStarts, SampleNames);



r=corr(intersect_iMutAF(bothmethods03),mutAF(bothmethods03,poscontrolsample))



%% Sputum 14 v Sputum 15 (Same patient seperated by 2 weeks)

good_at_either_time=(goodpos(:,2) | goodpos(:,3) | fixedpos(:,2) | fixedpos(:,3));
div_clickable_scatter_samplevsample(mutAF(:,2),  mutAF(:,3), 'Day 0 Mutation Allele Frequency', 'Day 0 Mutation Allele Frequency', 2, 3, good_at_either_time, strict_parameters, coveragethresholds, counts, fwindows, isolate_windows, cwindows, positions, mutations, RefGenome, ScafNames, ChrStarts, SampleNames);








%% Find genes mutated more than once across pateints and within samples

  %%Also find perecent overlap Within vs between samples -- Conclusion, not statistically siginificant

genelengththreshold=2000;
genelengththresholdmutator=1000;
  
%per sample
percentsame=zeros(NSample,1);
sample_genes=cell(7,1);
sample_genes_unique=cell(7,1);
sample_genes_fixed=cell(7,1);
sample_genes_fixed_unique=cell(7,1);
multiple_locations_within_sample=cell(7,1);
multiple_locations_within_sample_all=[];
multiple_locations_within_sampleperlength=cell(7,1);


for i=2:7
    
    
    %find diverse genes
    g=genes(goodpos(:,i) > 0 & (cellfun(@numel,genes) > 0)')'; gm=cell2mat(g); % gene names
    sample_genes{i}=str2num(gm(:,6:end));
    sample_genes_unique{i}=unique(sample_genes{i});
    
    %find fixed genes
    gm=cell2mat(genes(fixedpos(:,i)>0 & (cellfun(@numel,genes) > 0)')');
    sample_genes_fixed{i}=str2num(gm(:,6:end));
    sample_genes_fixed_unique{i}=unique(sample_genes_fixed{i});
    sample_genes_unique{i}(ismember(sample_genes_unique{i},sample_genes_fixed_unique{i}));
    
    [percentsame(i),multiple_locations_within_sample{i}] =div_pairwise_nts_in_same_gene([sample_genes{i}]); %;sample_genes_fixed{i}] );
    
    
    [percentsame(i),multiple_locations_within_sample{i}] =div_pairwise_nts_in_same_gene([sample_genes{i}]); %;sample_genes_fixed{i}] );
    
    sg=sample_genes{i};
    uniquegenes=unique(sg);
    mutspergene=zeros(ones,numel(uniquegenes));
    for k=1:numel(sg)
        j=find(uniquegenes==sg(k));
        mutspergene(j)=mutspergene(j)+1;
    end
    mutspergenelength=mutspergene./genesizes(uniquegenes);
    

    if mutator~=i
        multiple_locations_within_sampleperlength{i}=uniquegenes(mutspergenelength>(1/genelengththreshold) & mutspergene>1);
    else
        multiple_locations_within_sampleperlength{i}=uniquegenes(mutspergenelength>(1/genelengththresholdmutator) & mutspergene>1);
    end
    
    %only add nonmutators to "all" list -- this two tiered system allowed
    %inclusion of previous sample 
    
    if ismember(i,nonmutatorsamples)
        multiple_locations_within_sample_all=unique([multiple_locations_within_sample_all; multiple_locations_within_sampleperlength{i}]);
    end
        %  else
    %    multiple_locations_within_sample_all=unique([multiple_locations_within_sample_all; multiple_locations_within_sample{i}]);
   % end

end


%all loci
g=genes(sum(goodpos(:,3:end)+fixedpos(:,3:end),2) > 0  & (cellfun(@numel,genes) > 0)')'; gm=cell2mat(g); %all gene names
allg=str2num(gm(:,6:end));
percentsame_all=div_pairwise_nts_in_same_gene(allg);

numt=1000;
N=zeros(1000,1);
%bootstrap on all loci
for t=1:1000
     r=floor(rand(82,1).*numel(allg))+1;
     N(t)=div_pairwise_nts_in_same_gene(allg(r));
end
N=sort(N); disp([N(25) N(975)]);


%2 samples from same patient
g=genes(sum(goodpos(:,2:3)+fixedpos(:,2:3),2) > 0  & (cellfun(@numel,genes) > 0)')'; gm=cell2mat(g); %all gene names
percentsame_2sample=div_pairwise_nts_in_same_gene(str2num(gm(:,6:end)));



%for each patient, call a gene once. whats the overlap between samples?
ugenes=[];
for i=nonmutatorsamples
    ugenes=[ugenes; unique([sample_genes_unique{i}; sample_genes_fixed_unique{i}])];
end
%ugenes=[ugenes; unique([sample_genes_unique{2}; sample_genes_fixed_unique{2}; sample_genes_unique{3}; sample_genes_fixed_unique{3}])];
[gene_overlap_between_patients, multiple] =div_pairwise_nts_in_same_gene(ugenes);



%% Determine coverable positions and coverage --- perhaps faster if includes in creation of positions or counts, but not done currently

switch processed
    case 0
        savename=['callable_positions_' postfix];
        callablePos= div_genomic_positions_with_potential_to_call_diversity(filen,SampleDirs, strict_parameters, coveragethresholds, GenomeLength);
        %remove from callable if diverse in control
        callablePos(p(maf(:,1)<strict_parameters.min_control_MAF),:)=0;
        save(savename,'callablePos')
    case 1
        savename=['callable_positions_' postfix];
        load(savename)
end



%coverage positions in isolates generated with a python script and stored
%in a text file (based on FQ in strain.vcf)

 


%% Is it statistically unlikely to have this many genes with multiple mutations in a sample?

%genelengththreshold=2000;

numtrials=1000;

expectedNumberMultipleMutatedGenes=zeros(numtrials,size(counts,3));


for i=samples
    disp([i sum(goodpos(:,i))])
    for t=1:numtrials
        
        %get random positions on genome, restricting to those positions
        %which were callable
        f=find(callablePos(:,i)>0);
        genenums = bdolosa_genenums(f(floor(numel(f)*rand(sum(goodpos(:,i)),1))+1), cds, ChrStarts);
       
        sim_genes=unique(genenums);
        sim_mutspergene=zeros(1,numel(sim_genes));
        for k=1:numel(genenums)
            j=find(sim_genes==genenums(k));
            sim_mutspergene(j)=sim_mutspergene(j)+1;
        end
       
         sim_mutspergenelength=sim_mutspergene./genesizes(sim_genes);

         %this is a treshold criteria
     %   expectedNumberMultipleMutatedGenes(t,i)=sum(genesizes(sim_genes)<4000 & sim_mutspergene>1);

       expectedNumberMultipleMutatedGenes(t,i)=sum(sim_mutspergenelength>(1/genelengththreshold) & sim_mutspergene>1);        
    end
end

%Plot
figure(21); clf;
k=1;
for i=samples(1:end-1);
    
    subplot(numel(samples),1,k); hold on;

   
    [freq, bins]=hist(expectedNumberMultipleMutatedGenes(:,i),max(expectedNumberMultipleMutatedGenes(:,i)));
    bar(bins-.5,freq,'FaceColor', rgb('gray'), 'EdgeColor', 'w')
    
    bar(numel(multiple_locations_within_sampleperlength{i}),numtrials,.25,'FaceColor',rgb('DarkMagenta'),'EdgeColor', rgb('DarkMagenta'))
    
    %title(SampleNames(i).Sample)
      axis([-0.6 10 -.01 numtrials])
     set(gca,'XTick', 0:2:10)
    set(gca,'YTick', [.25*numtrials .5*numtrials .75*numtrials numtrials])
    set(gca,'YTicklabel', {'.25', '.50', '.75', '1'})
    k=k+1;
end;


%Plot hypermutator with slightly different settings
i=samples(end);
subplot(numel(samples),1,k); cla; hold on;
[freq, bins]=hist(expectedNumberMultipleMutatedGenes(:,i),max(expectedNumberMultipleMutatedGenes(:,i)));
bar(bins-.5,freq,1,'FaceColor', rgb('gray'), 'EdgeColor', 'w')
bar(numel(multiple_locations_within_sampleperlength{i}),numtrials,.25, 'FaceColor',rgb('DarkMagenta'),'EdgeColor', rgb('DarkMagenta'))
%title(SampleNames(i).Sample)
axis([-.6 10 -.001 numtrials/4])
set(gca,'YTick', [.05*numtrials .10*numtrials .15*numtrials .2*numtrials])
set(gca,'YTicklabel', {'.05', '.10', '.15', '.20'})
%set(gca,'XTick', 0:5:20)



observedngenes=cellfun(@numel,multiple_locations_within_sampleperlength)';
observedngenes=repmat(observedngenes,numtrials,1);


% 
% 
% 
% figure(749); clf; hold on; hist(othersAF,edges2);
% h=findobj(gca,'Type','patch'); set(h,'FaceColor', 'r', 'EdgeColor', 'w')
% xlabel('Major allele frequency'); ylabel('Number of positions passing strict quality filters')
% hist(selectionAF,edges2);
% legend({'Mutations in other genes', 'Mutation in genes under selection'}, 'Location', 'NorthWest')
% axis([0 1 -inf inf])
% %[a,b]=kstest2(selectionAF,othersAF);





%% Different way of asking the question-- expected and actual distribution of mutations per gene length


% this is commented out because its not simulation based


% %for now, just for patient 1
% figure(22); clf;
% k=0;
% for i=samples
%     
%     k=k+1;
%     observedmutatedgenes=gn((goodpos(:,i)>0 & gn>0));
%     observedmutatedgenesU=unique(observedmutatedgenes);
%     
%     %expected number of mutations for a gene is just that gene length divided
%     %by the size of the genome
%     [ehits,ebins]=hist([1000* numel(observedmutatedgenes) *genesizes(genesizes >0)/GenomeLength],50);
%     
%     
%     
%     %add each time a gene was found
%     observedmutspergene=zeros(1,numel(observedmutatedgenesU));
%     for i=1:numel(observedmutatedgenes)
%         j=find(observedmutatedgenesU==observedmutatedgenes(i));
%         observedmutspergene(j)=observedmutspergene(j)+1;
%     end
% 
%     
%     disp(max(observedmutspergene))
%    observedmutspergene=observedmutspergene./genesizes(observedmutatedgenesU);
%      [hits,bins]=hist([1 observedmutspergene],5000);   
% 
%   
%      figure(865); clf; hold on;
%    % bar(bins,hits/sum(hits),'r')
%     bar(ebins,ehits/sum(ehits),'b')
%     legend({'Observed genes', 'Expected - all genes'})
%     xlabel('Mutations per 1kb')
%     ylabel('Percent of genes')
% 
%     
%   %  axis([0 15 -inf inf])
%     
% %     subplot(numel(samples),1,k); hold on;
% % 
% %     hist([genesizes(observedmutatedgenes) 20000], 200)
% %     axis([0 14000 -inf inf])
% 
%     pause
% end


figure; hist(genesizes(genesizes>0),140)

%% Create expected dNdS

 


switch 1
    case 0
        
        %this section was built with a different strategy in mind.
        %cds_possibilites is the most important thing, and is transformed
        %below into the ever-useful percentN_types
        
        
        savename=['Expected_' postfix];
        
        [genomewide_possibilities, cds_possibilities] = div_mutation_type_probability_matrix(cds, GenomeLength, ChrStarts, sequences);
        [upper95dNdS, lower95dNdS, dNdStrialsmatrix, percentN_matrix, ~, ~, expectedI, expectedP, coding_strand_possibilities] = div_permute_expected(2000, combined_annotation_goodpos, genomewide_possibilities, cds_possibilities);
        save(savename, 'genomewide_possibilities', 'cds_possibilities', 'expectedI', 'expectedP', 'upper95dNdS', 'lower95dNdS', 'dNdStrialsmatrix', 'percentN_matrix') %, 'expected_intergenic', 'expected_promoter')
    case 1
        savename=['Expected_' postfix];
        load(savename)
end





%% Mutation-mutation probability matrix for each patient, and 6 classes of mutations

GCratio=2;


psmatrix=zeros(4,4,NSample); %polymorphism or substitution 
pImatrix=zeros(4,4,NSample); %only intragenic polymorphism 
sImatrix=zeros(4,4,1); %only intragenic substitutions -- count each only once!
imatrix=zeros(4,4,NSample); %polymorphism or substitution , only intragenic
mut2matrix=zeros(4,4,NSample); %polymorphism in multiple_locations_within_sample
omatrix=zeros(4,4,NSample); %polymorphism NOT in multiple_locations_within_sample, intragenic
pmatrix=zeros(4,4,NSample); %polymorphism, intragenic or intergenic


%written terribly cause substitutions used to be counted per patient-- now
%counted only once


for i=nonmutatorsamples
 
    pz=find(goodpos(:,i)+fixedpos(:,i)>0);
    for f=1:numel(pz)
        
        j=pz(f);
        original_index = ancnt(j);
        %find new nt
        if goodpos(j,i) > 0
            new_index = [maNT(j, i) minorNT(j, i)];
            new_index(new_index==original_index)=[]; %remove original
        else
            new_index = [ancnt(j) maNT(j, i)];
            
        end
        new_index(new_index==original_index)=[]; %remove original
        
        
        psmatrix(original_index, new_index,i)=psmatrix(original_index, new_index,i)+1;
        
        
        
        if (types(j)=='N' | types(j)=='S') & goodpos(j,i) > 0
            imatrix(original_index, new_index,i)=imatrix(original_index, new_index,i)+1;
            pImatrix(original_index, new_index,i)=pImatrix(original_index, new_index,i)+1;

            
            if find(multiple_locations_within_sample{i}==str2num(genes{j}(6:end))) 
                mut2matrix(original_index, new_index,i)=mut2matrix(original_index, new_index,i)+1;
            else
                omatrix(original_index, new_index,i)=omatrix(original_index, new_index,i)+1;
            end
        end
        
    end
    
    %AT
    %AC
    %AG transition
    %GC
    %GT
    %GA transition
    
end



%count intragenic and intergenic for all, including hypermutator
for i=samples
 
    pz=find(goodpos(:,i));
    for f=1:numel(pz)
        
        j=pz(f);
        original_index = ancnt(j);
        %find new nt
        if goodpos(j,i) > 0
            new_index = [maNT(j, i) minorNT(j, i)];
            new_index(new_index==original_index)=[]; %remove original
        else
            new_index = [ancnt(j) maNT(j, i)];
            
        end
        
        new_index(new_index==original_index)=[]; %remove original
        pmatrix(original_index, new_index,i)=pmatrix(original_index, new_index,i)+1;
    end
end

%positions where at least one patient has fixed mutation
pz=find(sum(fixedpos(:,nonmutatorsamples),2)>0); 
for f=1:numel(pz)
    j=pz(f);
    original_index = ancnt(j);
    new_index = unique(maNT(j,:));
    new_index(new_index==original_index)=[]; %remove original
    if types(j)=='N' | types(j)=='S'
        sImatrix(original_index, new_index)=sImatrix(original_index, new_index)+1;
    end
    
    
end




onlyImutO = div_matrix2_6types(imatrix);
onlySImutO = div_matrix2_6types(sImatrix);
onlyPImutO = div_matrix2_6types(pImatrix);
mut2mutO = div_matrix2_6types(mut2matrix);
mutO = div_matrix2_6types(psmatrix);
othersmutO= div_matrix2_6types(omatrix);
onlyPmutO = div_matrix2_6types(pmatrix);

percentN_types=div_matrix2_6types(percentN_matrix)./2;



figure(3002); clf; hold on;
NmutO=onlyPmutO; %NmutO(4:6,:)=NmutO(4:6,:) %./GCratio;
bar(1:numel(samples),NmutO(:,samples)');
legend({'AT->TA', 'AT->CG', 'AT->GC', 'GC->CG', 'GC->TA', 'GC->AT'})
ylabel('Number of SNPs')
set(gca,'Xtick',1:6)
set(gca,'Xticklabel',{'P1', 'P2', 'P3', 'P4', 'P5'})


figure(3005); clf; hold on;
NmutO=onlyPmutO; NmutO(4:6,:)=NmutO(4:6,:) %./GCratio;
bar(1:6,NmutO([1 2 4 5 3 6],samples));
%legend({'P1', 'P2', 'P3', 'P4', 'P5'})
ylabel('Number of SNPs')
set(gca,'Xtick',1:6)
set(gca,'Xticklabel',{'AT->TA', 'AT->CG', 'GC->CG', 'GC->TA', 'AT->GC', 'GC->AT'})
colormap([.85 .85 .85; .7 .7 .7; .5 .5 .5; .3 .3 .3; 0 0 0;])
%axis([0 7 0 190])


figure(3003); clf; hold on;
bar([1 2], [percentN_types'; 0 0 0 0 0 0], 'FaceColor', rgb('Gray'));
ylabel('Probability nonsynonymous')
set(gca,'Xtick',.67:.133:1.4)
set(gca,'Xticklabel',{'A->T', 'A->C', 'A->G', 'G->C', 'G->T', 'G->A'})
axis([0.5 1.5 0 1])


figure(3004); clf; hold on;
bar([1 2], [percentN_types'./(1-percentN_types'); 0 0 0 0 0 0], 'FaceColor', rgb('Gray'));
ylabel('Expected N/S')
set(gca,'Xtick',.67:.133:1.4)
set(gca,'Ytick',0:2:16)
set(gca,'Xticklabel',{'A->T', 'A->C', 'A->G', 'G->C', 'G->T', 'G->A'})
axis([0.5 1.5 0 15])




expected_percentN = (sum(onlyImutO(:,nonmutatorsamples),2)'*percentN_types) / sum(sum(onlyImutO(:,nonmutatorsamples),2)); 
expected = expected_percentN/(1 - expected_percentN); 


%% Types of mutations for each patient

figure(3004); clf; hold on;

sumtransversions=sum(onlyPmutO([1 2 4 5],:));
sumtransitions=sum(onlyPmutO([3 6],:));
h=bar([sumtransversions(samples); sumtransitions(samples)]','grouped')
%set(gca,'Xticklabel',{SampleNames(samples).Sample})
colormap([rgb('Grey'); rgb('Black')])
%set(h,'EdgeColor', 'k')
legend({'Transversions','Transitions'})
grid(gca)
set(gca,'Xtick',1:5)
set(gca,'Xticklabel',{'P1', 'P2', 'P3', 'P4', 'P5'})
ylabel('Polymorphic positions')
set(gca,'Ytick',0:50:350)
set(gca, 'GridLineStyle', ':');




%% transitions versus transversions -- only count in samples

notsamples=ones(size(goodpos,2),1);
notsamples(samples)=0;

%first get the types of mutations-- transformming matrix to an array for
%now
ntlist=zeros(numel(goodpos),2);
ntlist(goodpos>0,1)=maNT(goodpos>0);
ntlist(goodpos>0,2)=minorNT(goodpos>0);
ntlist(fixedpos>0,1)=maNT(fixedpos>0);
ntlist(fixedpos>0,2)=ancntm(fixedpos>0);

mutationtype=reshape(isTrTv(ntlist(:,1), ntlist(:,2)),size(goodpos,1),size(goodpos,2));
mutationtype(:,notsamples>0)=0;


%% Are transitions and transversions at different frequencies?

%major allele freq if major allele is mutation or fixed mutation,
%otherwise minor allle freq
TvAF=mutAF(mutationtype==1 & (goodpos>0 | fixedpos>0));
TrAF=mutAF(mutationtype==2 & (goodpos>0 | fixedpos>0));

[~,trtvafp]=kstest2(TrAF,TvAF);
disp(['Probability that transitions and transversion frequencies are drawn from same distribution: ' num2str(trtvafp)])


%In only hypermutator?
TvAFHM=mutAF(mutationtype(:,4)==1 & (goodpos(:,4)>0 | fixedpos(:,4)>0));
TrAFHM=mutAF(mutationtype(:,4)==2 & (goodpos(:,4)>0 | fixedpos(:,4)>0));

[~,trtvafHMp]=kstest2(TrAFHM,TvAFHM);
disp(['Probability that transitions and transversion frequencies are drawn from same distribution in hypermutator: ' num2str(trtvafHMp)])

%no hypermutator
amutationtype=mutationtype;
amutationtype(:,4)=0;
TvAFnHM=mutAF(amutationtype==1 & (goodpos>0 | fixedpos>0));
TrAFnHM=mutAF(amutationtype==2 & (goodpos>0 | fixedpos>0));
[~,trtvafnHMp]=kstest2(TvAFnHM,TrAFnHM);
disp(['Probability that transitions and transversion frequencies are drawn from same distribution in all except: ' num2str(trtvafnHMp)])



figure(561); clf; hold on;
hist(TrAF,50)
h=findobj(gca,'Type','patch');
set(h,'FaceColor', 'r', 'EdgeColor', 'w')
hist(TvAF,50)
legend({'Transitions' 'Transversions'})
xlabel('Frequency')
ylabel('Mutations in all samples')

figure(562); clf; hold on;
hist(TrAFHM,50)
h=findobj(gca,'Type','patch');
set(h,'FaceColor', 'r', 'EdgeColor', 'w')
hist(TvAFHM,50)
legend({'Transitions' 'Transversions'})
xlabel('Frequency')
ylabel('Mutations in mutator')


figure(563); clf; hold on;
hist(TvAFnHM,50)
h=findobj(gca,'Type','patch');
set(h,'FaceColor', 'r', 'EdgeColor', 'w')
hist(TrAFnHM,50)
legend({'Transversions','Transitions'})
xlabel('Frequency')
ylabel('Mutations not in hypermutator')


%[~,blah1]=kstest2(TrAFnHM,TrAFHM)
%[~,blah2]=kstest2(TvAFnHM,TvAFHM)



%% Are N and S at different frequencies?

typesmatrix(:,notsamples>0)=0;
atypesmatrix=typesmatrix;
atypesmatrix(:,4)=0;


NAF=mutAF(typesmatrix=='N' & (goodpos>0 | fixedpos>0));
SAF=mutAF(typesmatrix=='S' & (goodpos>0 | fixedpos>0));


step2=.03;
edges2=.03:step2:1;

figure(745); clf; hold on; hist(NAF,edges2);
h=findobj(gca,'Type','patch'); set(h,'FaceColor', 'r', 'EdgeColor', 'w')
xlabel('Major allele frequency'); ylabel('Number of positions passing strict quality filters')
hist(SAF,edges2);
legend({'N', 'S'}, 'Location', 'NorthWest')
axis([0 1 -inf inf])

[~,pns]=kstest2(SAF, NAF)

%% Are INTRAGENIC and INTERGENIC at difference frequencies?


cAF=mutAF((typesmatrix=='N'| typesmatrix=='S')  & (goodpos>0 | fixedpos>0));
ncAF=mutAF((typesmatrix=='I'| typesmatrix=='P') & (goodpos>0 | fixedpos>0));


step2=.03;
edges2=.03:step2:1;

figure(746); clf; hold on; hist(cAF,edges2);
h=findobj(gca,'Type','patch'); set(h,'FaceColor', 'r', 'EdgeColor', 'w')
xlabel('Major allele frequency'); ylabel('Number of positions passing strict quality filters')
hist(ncAF,edges2);
legend({'Coding', 'Noncoding'}, 'Location', 'NorthWest')
axis([0 1 -inf inf])

[~,codingp]=kstest2(cAF, ncAF)


%% Are genes under selection at difference frequencies?

othersAF=[];
selectionAF=[];

for i=nonmutatorsamples
    othersAF= [othersAF; mutAF(~ismember(gn,multiple_locations_within_sample{i}) & goodpos(:,i)>0,i)];
    selectionAF = [selectionAF; mutAF(ismember(gn,multiple_locations_within_sample{i}) & goodpos(:,i)>0,i)];
end

step2=.03;
edges2=.03:step2:1;

figure(749); clf; hold on; hist(othersAF,edges2);
h=findobj(gca,'Type','patch'); set(h,'FaceColor', 'r', 'EdgeColor', 'w')
xlabel('Major allele frequency'); ylabel('Number of positions passing strict quality filters')
hist(selectionAF,edges2);
legend({'Mutations in other genes', 'Mutation in genes under selection'}, 'Location', 'NorthWest')
axis([0 1 -inf inf])

[~,codingp]=kstest2(othersAF, selectionAF)


%% Intragenic versus promoter


nonpromAF=mutAF(typesmatrix=='I' & (goodpos>0 | fixedpos>0));
promAF=mutAF(typesmatrix=='P' & (goodpos>0 | fixedpos>0));


step2=.03;
edges2=.03:step2:1;

figure(747); clf; hold on; hist(nonpromAF,edges2);
h=findobj(gca,'Type','patch'); set(h,'FaceColor', 'r', 'EdgeColor', 'w')
xlabel('Major allele frequency'); ylabel('Number of positions passing strict quality filters')
hist(promAF,edges2);
legend({'Non-promoter intragenic', 'Promter intragenic'}, 'Location', 'NorthWest')
axis([0 1 -inf inf])

[~,pns]=kstest2(nonpromAF, promAF)


%% Intragenic v intergenic abundances


expectedRationIntragenic = sum(sum(sum(genomewide_possibilities(:,:,1:2))))/sum(sum(sum(genomewide_possibilities(:,:,1:4))));
observedRatioIntragenic = sum(ismember(typesmatrix(:),'NS')& goodpos(:))/sum(ismember(typesmatrix(:),'NSIP')& goodpos(:));

%% dN dS by transitions transversions

atypesmatrix=typesmatrix;
atypesmatrix(:,4)=0;
TvN = sum(typesmatrix=='N' & (goodpos>0 | fixedpos>0) & mutationtype==1);
TvS = sum(typesmatrix=='S' & (goodpos>0 | fixedpos>0) & mutationtype==1);
TrN = sum(typesmatrix=='N' & (goodpos>0 | fixedpos>0) & mutationtype==2);
TrS = sum(typesmatrix=='S' & (goodpos>0 | fixedpos>0) & mutationtype==2);


expected_Tv=(sum(percentN_types([1 2 4 5]))/4)/(1-(sum(percentN_types([1 2 4 5]))/4));
expected_Tr=(sum(percentN_types([3 6]))/2)/(1-(sum(percentN_types([3 6]))/2));






%% N S substitution(fixed) vs polymorphic (diverse)

pN = sum(sum(typesmatrix(:,nonmutatorsamples)=='N' & goodpos(:,nonmutatorsamples)>0));
pS = sum(sum(typesmatrix(:,nonmutatorsamples)=='S' & goodpos(:,nonmutatorsamples)>0));
sN = sum(types'=='N' & sum(fixedpos(:,nonmutatorsamples),2)>0);
sS = sum(types'=='S' & sum(fixedpos(:,nonmutatorsamples),2)>0);


mcdonaldkreitman_ish=fexact([pN pS; sN sS])


%% N S in multiple within samples

m2N=0;
m2S=0;
oN=0;
oS=0;
for i=nonmutatorsamples
    
    %only count mutations in genes evolved in parallel within patient
%     m2N= m2N+sum(ismember(gn,multiple_locations_within_sampleperlength{i}) & (types=='N')' & goodpos(:,i)>0);
%     m2S=m2S+sum(ismember(gn,multiple_locations_within_sampleperlength{i}) & (types=='S')' & goodpos(:,i)>0);
%     oN=oN+sum(~ismember(gn,multiple_locations_within_sampleperlength{i}) & (types=='N')' & goodpos(:,i)>0);
%     oS=oS+sum(~ismember(gn,multiple_locations_within_sampleperlength{i}) & (types=='S')' & goodpos(:,i)>0);
%     
    
    %do the same thing, but include mutations found in other patients in genes
    %under parallel evolution in another patient
    m2N= m2N+sum(ismember(gn,multiple_locations_within_sample_all) & (types=='N')' & goodpos(:,i)>0);
    m2S=m2S+sum(ismember(gn,multiple_locations_within_sample_all) & (types=='S')' & goodpos(:,i)>0);
    oN=oN+sum(~ismember(gn,multiple_locations_within_sample_all) & (types=='N')' & goodpos(:,i)>0);
    oS=oS+sum(~ismember(gn,multiple_locations_within_sample_all) & (types=='S')' & goodpos(:,i)>0);
    
end





%% expectations of N and S in fixed and diverse

expectedByPt=zeros(size(goodpos,2),1);
expectedByPtp=zeros(size(goodpos,2),1); %polymorphism
expectedByPts=zeros(size(goodpos,2),1); %substitution
expectedByPtm2=zeros(size(goodpos,2),1); %mutated twice within a sample

Nexpected=0;
Sexpected=0;
pNexpected=0;
pSexpected=0;
m2Nexpected=0;
m2Sexpected=0;
oNexpected=0;
oSexpected=0;

for i=nonmutatorsamples
    %either
    rN=sum(percentN_types.* onlyImutO(:,i));
    rS=sum(onlyImutO(:,i))-rN;
    Nexpected=Nexpected+rN;
    Sexpected=Sexpected+rS;
    expectedByPt(i)=rN/rS;
    
    %polymorphism
    rN=sum(percentN_types.* onlyPImutO(:,i));
    rS=sum(onlyPImutO(:,i))-rN;
    pNexpected=pNexpected+rN;
    pSexpected=pSexpected+rS;
    expectedByPtp(i)=rN/rS;
    
    %multiple
    rN=sum(percentN_types.* mut2mutO(:,i));
    rS=sum(mut2mutO(:,i))-rN;
    m2Nexpected=m2Nexpected+rN;
    m2Sexpected=m2Sexpected+rS;    
    
    %notmultiple
    rN=sum(percentN_types.* othersmutO(:,i));
    rS=sum(othersmutO(:,i))-rN;
    oNexpected=oNexpected+rN;
    oSexpected=oSexpected+rS;    
    
    
       
end


%substitution -- count once
sNexpected=sum(percentN_types.* onlySImutO);
sSexpected=sum(onlySImutO)-sNexpected;


expectedp=pNexpected/pSexpected;
expecteds=sNexpected/sSexpected;
expectedm2=m2Nexpected/m2Sexpected;
expectedothers=oNexpected/oSexpected;





[upolymorphic,lpolymorphic]=binomialCIdNdS(pN,pS,expectedp);
[usubstitution,lsubstitution]=binomialCIdNdS(sN,sS,expecteds);
[um2,lm2]=binomialCIdNdS(m2N,m2S,expectedm2);
[uothers,lothers]=binomialCIdNdS(oN,oS,expectedothers);


figure(997); clf; hold on;

bar(0, pN/pS/expectedp, 'FaceColor', tamismagenta, 'EdgeColor', 'none', 'BarWidth', .65);
plot([0 0],[lpolymorphic upolymorphic], '-', 'Color', rgb('DarkGray'), 'LineWidth',1);

bar(1, sN/sS/expecteds, 'FaceColor', tamisgreen, 'EdgeColor', 'none','BarWidth', .65); 
plot([1 1],[lsubstitution  usubstitution],  '-', 'Color', rgb('DarkGray'), 'LineWidth',1);

plot([-1 5], [1 1], 'k-', 'LineWidth',1) %expected


axis([-.7 1.5 0.65 5 ])
set(gca,'YScale','log')
set(gca,'YTick', [.5 1 2 3 4 8 16])
%set(gca,'YTicklabel', {'.7', '.8', '1', '1.2', '1.4' '1.7', '2'})
ylabel('dN/dS (log scale)')
set(gca,'Xtick',[0 1])
set(gca,'Xticklabel',{'Polymorphisms','Substitutions'} )



figure(998); clf; hold on;

bar(0, m2N/m2S/expectedm2, 'FaceColor', rgb('Red'), 'EdgeColor', 'none', 'BarWidth', .65);
plot([0 0],[lm2 um2], '-', 'Color', rgb('DarkGray'), 'LineWidth',1);
plot([-1 5], [1 1], 'k-', 'LineWidth',1) %expected


bar(1, oN/oS/expectedothers, 'FaceColor', rgb('Black'), 'EdgeColor', 'none', 'BarWidth', .65);
plot([1 1],[uothers lothers], '-', 'Color', rgb('DarkGray'), 'LineWidth',1);
plot([-1 5], [1 1], 'k-', 'LineWidth',1) %expected


axis([-.7 1.5 0.65 6 ])
set(gca,'YScale','log')
set(gca,'YTick', [1 2 3 4 5])
ylabel('dN/dS (log scale)')
set(gca,'Xtick',[0 1])
%set(gca,'Xticklabel',{'Genes mutated multiple times within the same sample','Others'} )




% same plot, hortizonal 
figure(998); clf; hold on;

barh(1, m2N/m2S/expectedm2, 'FaceColor', rgb('DarkMagenta'), 'EdgeColor', 'none', 'BarWidth', .65);
plot([lm2 um2], [1 1], '-', 'Color', rgb('DarkGray'), 'LineWidth',1);


barh(0, oN/oS/expectedothers, 'FaceColor', rgb('Black'), 'EdgeColor', 'none', 'BarWidth', .65);
plot([uothers lothers], [0 0], '-', 'Color', rgb('DarkGray'), 'LineWidth',1);

plot( [1 1], [-1 5], 'k-', 'LineWidth',1) %expected


axis([0.9 9 -.7 1.5 ])
set(gca,'XScale','log')
set(gca,'XTick', [1 2  4 6 8])
xlabel('dN/dS (log scale)')
set(gca,'Ytick',[0 1])
set(gca,'Yticklabel',{' ',''} )
%set(gca,'Xticklabel',{'Genes mutated multiple times within the same sample','Others'} )






 %% Random dNdS-- confirm dNdS works -- when have time rerun this and save the histograms
 %also good for p values
% 
% numtrials=1000;
% dNdSrandom=zeros(numtrials,1);
% dNdSrandom2=zeros(numtrials,1);
% 
% flatimatrix=sum(imatrix(:,:,nonmutatorsamples),3);
% nullmatrix=25*ones(4,4); flatimatrix(1,1)=0; flatimatrix(2,2)=0; flatimatrix(3,3)=0; flatimatrix(4,4)=0;
% 
% 
% 
% for i=1:numtrials
% 
%     %generate mutations of same distribution as imatrix randomly on genome
%     
%     randp=floor(GenomeLength*rand(3*sum(flatimatrix(:)),1)); %grab three times as many positions because we need to pick positions with regard to imatrix
%     [~,  ~, randmutations, ~] = annotate_mutations_auto_gb(p2chrpos(randp,ChrStarts),ScafNames,RefGenome) ;
%     
%     
%     %calculate N and S
%     N=0;
%     S=0;
%     for old=1:4
%         k=0;
%         nt=NTs(old);
%         nt_mut_index=find([randmutations.ref]==nt & cellfun(@numel,{randmutations(:).AA})>0); %intragenic and ref is proper nt
%         randmutations_nt=randmutations(nt_mut_index(floor(numel(nt_mut_index)*rand(sum(flatimatrix(old,:)),1)+1)));
%         for new=1:4
%             if old~=new
%                 aas=cell2mat({randmutations_nt(k+1:k+flatimatrix(old,new)).AA}');
%                 N=N+sum(aas(:,old)~=aas(:,new));
%                 S=S+sum(aas(:,old)==aas(:,new));
%                 k=k+flatimatrix(old,new);
%             end
%         end
%     end
%     dNdSrandom(i)=N/S;
%   
%     
%     %reapeat for nullmatrix
%         
%     %calculate N and S
%     N=0;
%     S=0;
%     for old=1:4
%         k=0;
%         nt=NTs(old);
%         nt_mut_index=find([randmutations.ref]==nt & cellfun(@numel,{randmutations(:).AA})>0); %intragenic and ref is proper nt
%         randmutations_nt=randmutations(nt_mut_index(floor(numel(nt_mut_index)*rand(sum(nullmatrix(old,:)),1)+1)));
%         for new=1:4
%             if old~=new
%                 aas=cell2mat({randmutations_nt(k+1:k+nullmatrix(old,new)).AA}');
%                 N=N+sum(aas(:,old)~=aas(:,new));
%                 S=S+sum(aas(:,old)==aas(:,new));
%                 k=k+flatimatrix(old,new);
%             end
%         end
%     end
%     dNdSrandom2(i)=N/S;
%     
%     
% end
% 
%     
    
   




%% Plot MAF along genome


figure(311); clf; hold on;
set(311,'Position',[0 0 scrsz(3)/2 scrsz(3)]);clf;hold on;

approx_sizes=[3413073,2176129,831198];
colors='bgrkym';
subplotsizes=[1 3 10];

for j=1:3
    
    toplot=zeros(numel(p), NSample-1);
    x=sum(goodpos,2)>0 & positions(:,1)==j;
    for i=samples
        toplot(goodpos(:,i)>0 & positions(:,1)==j,i-1)= 1 - maf(goodpos(:,i)>0& positions(:,1)==j,i);
    end
    
    
    subplot(3,subplotsizes(j),(subplotsizes(j)*(j-1))+1:(subplotsizes(j)*(j-1))+j)
    
    bar(positions(x,2), toplot(x,:), 10, 'grouped')
    colormap jet
    
    title(['Chromosome ' num2str(j)])
    legend({SampleNames(2:end).Sample}, 'Position')
    
    % xlabel('Position on chromsome')
    ylabel('Minor allele frequency')
    axis([0 approx_sizes(j) 0 .5])
    
end


% % %% This piece of code only works for a different set of things re:
%false positives
% % 
% % 
% % figure(312); clf; hold on;
% % set(312,'Position',[0 0 scrsz(3)/2 scrsz(3)]);clf;hold on;
% % 
% % approx_sizes=[3413073,2176129,831198];
% % colors='kkkkkkk';
% % subplotsizes=[1 3 10];
% % 
% % for j=1:3
% % 
% %     subplot(3,subplotsizes(j),(subplotsizes(j)*(j-1))+1:(subplotsizes(j)*(j-1))+j)
% %     
% %     y=zeros(approx_sizes(j),1);
% %     y(positions(2,sum(minorAF(positions(1,:)==j,:)>.03,2)>0))=1;
% %     
% %     
% %     plot(smooth(y,1000),'.')
% %     colormap jet
% %     
% %     title(['Chromosome ' num2str(j)])
% %    % legend({SampleNames(2:end).Sample}, 'Position')
% %     
% %     % xlabel('Position on chromsome')
% %  %   ylabel('Minor allele frequency')
% %  %   axis([0 approx_sizes(j) 0 .5])
% %     
% % end
% % 



%% Histograms of MAF as a function of N, S, or I

step=.025;
edges=.5:step:1;
allN=zeros(numel(edges),1);
allS=zeros(numel(edges),1);
allI=zeros(numel(edges),1);
allP=zeros(numel(edges),1);

%each sample & all
for i=nonmutatorsamples
    figure(850+i);clf;
    N=reshape(histc(maf(goodpos(:,i)>0 & (types=='N')',i),edges),size(edges'));
    S=reshape(histc(maf(goodpos(:,i)>0 & (types=='S')',i),edges),size(edges'));
    I=reshape(histc(maf(goodpos(:,i)>0 & (types=='I')',i),edges),size(edges'));
    P=reshape(histc(maf(goodpos(:,i)>0 & (types=='P')',i),edges),size(edges'));
    bar(edges+step/2, [N S I P],'stacked');
    xlabel('Major allele frequency')
    ylabel('Number of positions passing strict quality filters')
    title(SampleNames(i).Sample)
    axis([.5 1 -inf inf])
    legend('Nonsynonymous', 'Synonymous','Intergenic','Promoter', 'Location', 'NorthWest')
    colormap cool
    allN=allN+N;
    allS=allS+S;
    allI=allI+I;
    allP=allP+P;
    
end



%all
figure(860);clf; hold on;
bar(edges+step/2, [allN allS allI allP],'stacked');
xlabel('Major allele frequency')
ylabel('Number of positions passing strict quality filters')
title('All nonmutatorsamples')
axis([.5 1 -inf inf])
legend('Nonsynonymous', 'Synonymous','Intergenic','Promoter', 'Location', 'NorthWest')
colormap cool

%isolates
figure(861);clf;
N=reshape(histc(iMAF(itypes=='N'),edges),size(edges'));
S=reshape(histc(iMAF(itypes=='S'),edges),size(edges'));
I=reshape(histc(iMAF(itypes=='I'),edges),size(edges'));
P=reshape(histc(iMAF(itypes=='P'),edges),size(edges'));
bar(edges+step/2, [N S I P],'stacked');
xlabel('Major allele frequency')
ylabel('Number of positions')
title('S2 isolates')
axis([.5 1 -inf inf])
legend('Nonsynonymous', 'Synonymous','Intergenic','Promoter', 'Location', 'NorthWest')
colormap cool


%% Are N and S in different places in genes?

figure(723); clf; hold on;
Ndist=[];
Sdist=[];
for i=nonmutatorsamples
    Ndist=[Ndist [combined_annotation_all(types(goodpos(:,i) > 0)=='N').nt_pos]];
    Sdist=[Sdist [combined_annotation_all(types(goodpos(:,i) > 0)=='S').nt_pos]];
end

hist(Ndist,50)
h=findobj(gca,'Type','patch');
set(h,'FaceColor', 'r', 'EdgeColor', 'w')
hist(Sdist,50)
legend({'N' 'S'})
xlabel('nt pos in gene')
ylabel('times found')









%% Fixed mutations

%Fails McDonald-Kreitman test just barely
% Is significantly positive dNdS if you do not include mutator sample

insamples=repmat(ismember(1:7,nonmutatorsamples),size(goodpos,1),1);



Nd=typesmatrix =='N' & goodpos>0 & insamples;
Sd=typesmatrix =='S' & goodpos>0 & insamples;

Nf=typesmatrix =='N' & fixedpos>0 & insamples;
Sf=typesmatrix =='S' & fixedpos>0 & insamples;


Nfixed=sum(Nf(:));
Sfixed=sum(Sf(:));
[uf,lf]=binomialCIdNdS(Nfixed,Sfixed,expected);


NAll=sum(Nd(:));
SAll=sum(Sd(:));


[uAll,lAll]=binomialCIdNdS(NAll,SAll,expected);

withinPt=sum(Nd)./sum(Sd)/expected;
betweenPt=sum(Nf)./(sum(Sf)+.01)/expected;

withinPtAll=sum(Nd(:))./sum(Sd(:))/expected;
[uw,lw]=binomialCIdNdS(sum(Nd(:)),sum(Sd(:)),expected);

betweenPtAll=sum(Nf(:))./sum(Sf(:))/expected;
[ubtwn,lbtwn]=binomialCIdNdS(sum(Nf(:)),sum(Sf(:)),expected);


a=fexact([sum(Nd(:)) sum(Sd(:)); sum(Nf(:)) sum(Sf(:))]);

[a,b]=kstest2(maf(Nd>0), maf(Sd>0));

for i=nonmutatorsamples
    within=sum(Nd(:,i))/sum(Sd(:,i))/expectedByPt(i);
    between=sum(Nf(:,i))/(sum(Sf(:,i))+.01)/expectedByPt(i);
    both=sum(sum(Nd(:,i))+sum(Nf(:,i)))/(sum(Sd(:,i))+sum(Sf(:,i)))/expectedByPt(i);
    
    [uboth,lboth]=binomialCIdNdS(sum(sum(Nd(:,i))+sum(Nf(:,i))),sum(sum(Sd(:,i))+sum(Sf(:,i))),expectedByPt(i));
    disp([within between both lboth uboth]);
    
   [~, kstestp]= kstest2(maf((types=='N' & (goodpos(:,i)>0)')>0,i)',maf((types=='S' & (goodpos(:,i)>0)')>0,i)');

  %  figure; clf; hold on; hist(maf((types=='N' & (goodpos(:,i)>0)')>0,i)'); hist(maf((types=='S' & (goodpos(:,i)>0)')>0,i)');
 %   a=fexact([sum(Nd(:,i)) sum(Sd(:,i)); sum(Nf(:,i)) sum(Sf(:,i))])

end




%% Frequencies of genes diverse in more than one patient

goodmaf=maf; goodmaf(goodpos<1)=nan;




%All goodpositions
figure(322); clf; hold on; hist(goodmaf(goodpos>0),50);
h=findobj(gca,'Type','patch'); set(h,'FaceColor', 'r', 'EdgeColor', 'w')
xlabel('Major allele frequency'); ylabel('Number of positions passing strict quality filters')
axis([.5 1 -inf inf])

%nt positions diverse/fixed in more than one patient
nts_in_multiple_pts=(sum(goodpos(:,3:end)+fixedpos(:,3:end),2)>1);
nts_in_multiple_pts_mat=repmat(nts_in_multiple_pts,1,7);
hist(goodmaf(nts_in_multiple_pts_mat),50);
legend({'All', 'NTs diverse in more than one patient'})

ntmaf=goodmaf(nts_in_multiple_pts_mat);
[~,gp]=kstest2(ntmaf(ntmaf>0), goodmaf(goodpos>0));

%nt positions in genes that are diverse in more than one patient
mgenes_maf=[];
p_mgenes_maf=[];
mgenes_ind=[];
p_mgenes_ind=[];
for i=1:numel(genes)
    if numel(genes{i} > 0) & ismember(str2num(genes{i}(6:end)),multiple)
        s=goodmaf(i,:); s=s(s>0);
        if size(s)>0
            mgenes_ind(end+1)=i;
            mgenes_maf(end+1:end+numel(s))=s;
        end
    end
    if numel(genes{i} > 0) & ismember(str2num(genes{i}(6:end)),multiple_locations_within_sample_all)
        s=goodmaf(i,:); s=s(s>0);
        if size(s)>0
            p_mgenes_ind(end+1)=i;
            p_mgenes_maf(end+1:end+numel(s))=s;
        end
    end
end


figure(323); clf; hold on; hist(goodmaf(goodpos>0),50);
h=findobj(gca,'Type','patch'); set(h,'FaceColor', 'r', 'EdgeColor', 'w')
xlabel('Major allele frequency'); ylabel('Number of positions passing strict quality filters')
axis([.5 1 -inf inf])
hist(mgenes_maf,50);
legend({'All', 'Genes diverse in more than one patient'})
[~,mgenesp]=kstest2(mgenes_maf, goodmaf(goodpos>0))


figure(324); clf; hold on; hist(goodmaf(goodpos>0),50);
h=findobj(gca,'Type','patch'); set(h,'FaceColor', 'r', 'EdgeColor', 'w')
xlabel('Major allele frequency'); ylabel('Number of positions passing strict quality filters')
axis([.5 1 -inf inf])
hist(p_mgenes_maf,50);
legend({'All', 'Genes with more than one position diverse within a patient'})
[~,p_mgenes_p]=kstest2(p_mgenes_maf', goodmaf(goodpos>0))





%% Regions with more than one diverse position in 50bp window in same individual



switch processed
    case 0
        savename=['NearbyNTpairs'];
        ntpairing= div_analyze_nearby_nts(p, goodpos, GenomeLength, SampleNames, SampleDirs, ChrStarts, combined_annotation_all, 50);
        save(savename, 'ntpairing')
    case 1
        savename=['NearbyNTpairs'];
        load(savename)
end

NTs='ATCG';

%display by position
figure(260);clf;
numpairs=numel(ntpairing);
wts=zeros(numpairs,1);
alt1s=zeros(numpairs,1);
alt2s=zeros(numpairs,1);
doubles=zeros(numpairs,1);

total=zeros(numpairs,1);
%expected=zeros(numpairs,1);


titles={};

for i=1:numpairs
    
    subplot(4,floor(numpairs/4)+1,i); hold on;
    
    m=ntpairing(i).hits; ref1=ntpairing(i).ref1; ref2=ntpairing(i).ref2;
    treads=sum(m(:));
    
    %ref, ref
    
    wt=m(ref1,ref2); %save and plot later to make legend easy to plot
    wts(i)=wt/treads; %save for later use
    %plot ref
    bar(1,wt,'k')
    
    
    %save items for calcuation of double mutant
    total(i)=sum(m(:));
    %expected(i)=(sum(m(:))-sum(m(ref1,:)))*(sum(m(:))-sum(m(ref2,:)))/(sum(m(:))*sum(m(:)));
    
    
    m(ref1,ref2)=0; blank=[0 0 0 0]; %remove for easy calculation

    [~,alt2]=sort(m(ref1,:));
    alt2=alt2(end);
    [~,alt1]=sort(m(:,ref2));
    alt1=alt1(end);
    
    %ref, alt
    bar([2 1],[m(ref1,:); blank],'stacked');
    alt1s(i)=sum(m(ref1,:))/treads;
    
    %alt, ref
    bar([3 2],[m(:,ref2)'; blank],'stacked');
    alt2s(i)=sum(m(:,ref2))/treads;
    

    
    
    %alt, alt
    m(ref1,:)=0; m(:,ref2)=0;
    bar(4,sum(m(:)),'g');
    doubles(i)=sum(m(:))/treads;
    
    % ylabel('Reads supporting')
    set(gca,'Xtick',1:4)
    set(gca,'Xticklabel',{'A', '1', '2', 'D'})


    if size(combined_annotation_all(ntpairing(i).index1).muts) > 0
        mut1=combined_annotation_all(ntpairing(i).index1).muts{1};
        mut2=combined_annotation_all(ntpairing(i).index2).muts{1};
    else
        mut1=[NTs(ref1) '->' NTs(alt1)];
        mut2=[NTs(ref2) '->' NTs(alt2)];
    end

    
    % set(gca,'Xticklabel',{[NTs(ref1) '-' NTs(ref2)], [NTs(ref1) '-X'], ['X-' NTs(ref2)], 'X-X'})
    % [NTs(ref1) ' - ' NTs(ref2) '  (wt)'], [NTs(ref1) ' - X'], ['X - ' NTs(ref2)], 'Double mut'} )
    
    %set(gca,'Xticklabel',{[NTs(ref1) ' - ' NTs(ref2) '  (wt)'], [NTs(ref1) ' - X'], ['X - ' NTs(ref2)], 'Double mut'} )
    
    c1=p2chrpos(ntpairing(i).pos1, ChrStarts);  c2=p2chrpos(ntpairing(i).pos2, ChrStarts);
    
    if isempty(strfind(ntpairing(i).annotation1(1),'B'))
        t1='Intergenic';
    else
        x=strfind(ntpairing(1).annotation1,' ');
        t1=ntpairing(i).annotation1(1:x(1));
    end
    title({[ntpairing(i).sample ': '   t1]})
    %title({[ntpairing(i).sample ':    Chr ' num2str(c1(1))] ;[ num2str(c1(2)) ': ' ntpairing(i).annotation1] ; [ num2str(c2(2)) ': '  ntpairing(i).annotation2] })
    titles{end+1}=[ntpairing(i).sample ': '   t1];
    
    
%     %display mathematica code
     disp(['g' num2str(i) '=mycylinders[' num2str(wts(i)) ',' num2str(alt1s(i)) ',' num2str(alt2s(i)) ',' num2str(doubles(i)) ', "' mut1 '","' mut2 '","'   titles{end} '"];' ])
end



%display by type
figure(261);clf; hold on;

h=barh([wts alt1s alt2s doubles], 'stacked')
colormap([1 1 1; .7 .1 .1; 0 .5 1; .7 .1 .7;])
set(gca,'Ytick',1:9)
set(gca,'Yticklabel',titles)


%model
figure(262);clf; hold on;

h=barh([.4 .3 .3 0; .7 0 0 .3; .3 .304 .304 .304^2], 'stacked')
colormap([1 1 1; .7 .1 .1; 0 .5 1; .7 .1 .7;])
set(gca,'Ytick',1:3)
%legend({'No mutation', 'Mutant 1', 'Mutant 2', 'Double Mutant'})
axis([0 1 0 4])








%% Average # of SNPs between isolates == 23



%% Uncategorized analyses

ntmutated=[];
ntdiverse=[];
tmutated=[];
smutated=[]; %how many samples is it mutated in?
pmutated=[]; %how many patients is it mutated in? to make this a fair comparision between patients, just check once
pdiverse=[];
pfixed=[];


for i=1:numel(multiple);
    ntmutated(i)=sum(gn==multiple(i) & sum(goodpos(:,samples)+fixedpos(:,samples),2)>0);
    ntdiverse(i)=sum(gn==multiple(i) & sum(goodpos(:,samples),2)>0);
    tmutated(i)=sum(sum(goodpos(gn==multiple(i),samples)+fixedpos(gn==multiple(i),samples),2));
    smutated(i)=sum(sum(goodpos(gn==multiple(i),samples)+fixedpos(gn==multiple(i),samples))>0);
    pmutated(i)=sum(sum(goodpos(gn==multiple(i),3:end)+fixedpos(gn==multiple(i),3:end))>0);
    pdiverse(i)=sum(sum(goodpos(gn==multiple(i),3:end))>0);
    pfixed(i)=sum(sum(fixedpos(gn==multiple(i),3:end))>0);

end




fgenes=unique(gn(sum(fixedpos(:,3:end),2)>0));
dgenes=unique(gn(sum(goodpos(:,3:end),2)>0));
mutated_genes=unique(gn(sum(fixedpos(:,3:end)+goodpos(:,3:end),2)>0));



nonanc=zeros(numel(multiple_locations_within_sample_all),numel(samples));
gnmat=repmat(gn,1,numel(samples));
for i=1:numel(multiple_locations_within_sample_all);
    
    nonanc(i)=sum(gn==mutated_genes(i) & sum(goodpos(:,samples)+fixedpos(:,samples),2)>0);
    
    ntmutated(i)=sum(gn==mutated_genes(i) & sum(goodpos(:,samples)+fixedpos(:,samples),2)>0);
    ntdiverse(i)=sum(gn==mutated_genes(i) & sum(goodpos(:,samples),2)>0);
    tmutated(i)=sum(sum(goodpos(gn==mutated_genes(i),samples)+fixedpos(gn==mutated_genes(i),samples),2));
    smutated(i)=sum(sum(goodpos(gn==mutated_genes(i),samples)+fixedpos(gn==mutated_genes(i),samples))>0);
    pmutated(i)=sum(sum(goodpos(gn==mutated_genes(i),3:end)+fixedpos(gn==mutated_genes(i),3:end))>0);
    pdiverse(i)=sum(sum(goodpos(gn==mutated_genes(i),3:end))>0);
    pfixed(i)=sum(sum(fixedpos(gn==mutated_genes(i),3:end))>0);
end

for i=1:numel(multiple);

end

specific_pressures=multiple(ntmutated==1);



% %% Find which genes to get fasta files for and which samples to get from
% 
% either=goodpos+fixedpos;
% 
% unknownp=find(types'=='U' & sum(either,2)>0);
% [unknown_genes,up]=unique(genes(unknownp));
% unknownp=unknownp(up);
% 
% unknown_strain=zeros(numel(unknownp),1);
% for i=1:numel(unknownp) 
%     %skip possibility of using isogenic control
%     unknown_strain(i) = find(either(unknownp(i),2:end)==0,1)+1;
% end
% 
% %save('unknown_gene_list','unknown_genes','unknownp','unknown_strain')


%% Summary statistics, how many positions met criteria?





%% Early vs late in isolates -- not significant

%For now, directly specify each branch in finding early mutations

NEarly=0;
SEarly=0;
anctntI=ancnt(ismember(allp,ip));
typesI=types(ismember(allp,ip));

%first two in triplet is ingroup, last is an outgroup
pairs={['TL033'; 'TL040'], ['TL030'; 'TL025'], ['TL035'; 'TL029']};
%['TL001'; 'TL013'], ['TL002'; 'TL004']}; 

for i=1:numel(pairs)
    ind1=find(strcmp(pairs{i}(1,:), {isolateNames.Sample}));
    ind2=find(strcmp(pairs{i}(2,:), {isolateNames.Sample}));
        
    uniquep = (isolateCalls(:,ind1)==isolateCalls(:,ind2)) & (isolateCalls(:,ind1)~=NTs(anctntI)'); 
    disp([ind1, ind2, sum(uniquep)])
    NEarly = NEarly+ sum((typesI =='N')' & uniquep==1);
    SEarly = SEarly+ sum((typesI =='S')' & uniquep==1);
end
% 
% %Unique to TL038
% ind1=find(strcmp('TL038', {isolateNames.Sample}));
% NEarly = NEarly+ sum((typesI =='N')' & (isolateCalls(:,ind1)~=NTs(anctntI)'));
% SEarly = SEarly+ sum((typesI =='S')' & (isolateCalls(:,ind1)~=NTs(anctntI)'));
% 


%% Plot mutational grid for 1st patient



treeorder=['TL005'; 'TL030'; 'TL028'; 'TL027'; 'TL018'; 'TL007'; 'TL021'; 'TL006'; 'TL014'; 'TL009'; 'TL008'; 'TL038'; 'TL016'; 'TL017'; 'TL010'; 'TL026'; 'TL015'; 'TL025'; 'TL004'; 'TL002'; 'TL029'; 'TL035'; 'TL041'; 'TL001'; 'TL013'; 'TL033'; 'TL040'; 'TL003'; 'TL037'];
[~,isolateorder]=ismember({isolateNames.Sample},treeorder);

figure(86)
clf; hold on;
numisolates=size(isolateCalls,2);
ancntImat=repmat(ancnt(ismember(allp,ip)),1,29);
typesImat=repmat(types(ismember(allp,ip)),29,1)';

%Plot N
[Nx,Ny]=find(isolateCalls~=NTs(ancntImat) & isolateCalls~='N' & typesImat=='N');
plot(Nx./2,1+numisolates-isolateorder(Ny),'r+','MarkerSize', 18)

%Plot S
[Sx,Sy]=find(isolateCalls~=NTs(ancntImat) & isolateCalls~='N' & typesImat=='S');
plot(Sx./2,1+numisolates-isolateorder(Sy),'b+','MarkerSize', 18)

%Plot I and P
[Ix,Iy]=find(isolateCalls~=NTs(ancntImat) & isolateCalls~='N' & typesImat~='N' & typesImat~='S');
plot(Ix./2, 1+numisolates-isolateorder(Iy),  'LineStyle', 'none', 'Marker', '+', 'MarkerEdgeColor', rgb('Gray'),'MarkerSize', 18)


%Plot lines
plot([-5*ones(1,numisolates); (size(isolateCalls,1)/2)*ones(1,numisolates)+5],[1:numisolates; 1:numisolates] ,'k','LineStyle','-')
axis([-3 95 -5 32])


%% Pairwise distances between strains within and between patients

%Within patients
diff=zeros(size(samples));
for i=samples
    ps=maf(goodpos(:,i)>0,i);
    qs=minorAF(goodpos(:,i)>0,i);
    diff(i)=sum(2*ps.*qs); 
end


avgdistancefromanc=zeros(1,NSample);
avgdistancefromPTanc=zeros(1,NSample);
numfixed=zeros(1,NSample);


for i=1:NSample
    nonancf=mutAF(:,i);
    avgdistancefromPTanc(i)=sum(nonancf(goodpos(:,i)>0));
    numfixed(i)=sum(fixedpos(:,i)>0 & ancnt~=maNT(:,i));
    avgdistancefromanc(i)=avgdistancefromPTanc(i)+ numfixed(i);
end


avgdistancefromothers=zeros(size(samples));
for i=samples
    others=samples;
    others(samples==i)=[];
    avgdistancefromothers(i)=mean(avgdistancefromanc(others))+mean(sum(fixedpos(:,others)))+avgdistancefromanc(i)+sum(fixedpos(:,i));
end


figure(421); clf; hold on;
bar([avgdistancefromPTanc(samples)./avgdistancefromanc(samples); numfixed(samples)./avgdistancefromanc(samples)]','stacked')
colormap([tamismagenta; tamisgreen;])
ylabel('Percent of SNPs per isolate relative to LCA of 5 patients')
%legend({'Within-patient Polymorphism', 'Substitution'})
set(gca,'Xtick',1:numel(samples))
set(gca,'Xticklabel',{'P1', 'P2', 'P3', 'P4', 'P5'})
axis([0 numel(samples)+1 0 1])


figure(483); clf; hold on;
h=bar(diff(samples), .5);
ylabel('Mean # SNPs between two isolates from patient')
set(gca,'Xtick',1:numel(samples))
set(gca,'Xticklabel',{SampleNames(samples).Sample})
set(h,'FaceColor', 'w', 'EdgeColor', rgb('DarkRed'), 'LineWidth', 3)


figure(484); clf; hold on;
h=barh([4 3 2 1], avgdistancefromPTanc(nonmutatorsamples)*GenomeLength./sum(callablePos(:,nonmutatorsamples)), .5);
disp(avgdistancefromPTanc(nonmutatorsamples)*GenomeLength./sum(callablePos(:,nonmutatorsamples))/2.1)
%xlabel('Average number of SNPs to patient LCA (?p)')
set(gca,'Ytick',1:numel(samples(1:end-1)))
set(gca,'Yticklabel',{'P4', 'P3', 'P2', 'P1'})
% set(gca,'XAxisLocation','Top')
% set(gca, 'Xtick', 0:5.25:21)
% set(gca, 'Xticklabel', [0:5.25:21]/2.1)
% xlabel('Average years to patient LCA (?p)')
set(h,'FaceColor', rgb('Black'), 'EdgeColor', 'None')
axis([0 28 0 numel(samples(1:end-1))+1 ])
grid(gca)
set(gca, 'GridLineStyle', ':');


figure(485); clf; hold on;
barh([sum(goodpos(:,samples))./sum(goodpos(:,samples)+fixedpos(:,samples));sum(fixedpos(:,samples))./sum(fixedpos(:,samples)+goodpos(:,samples))]','stacked')
xlabel('Fraction of mutations')
set(gca,'Ytick',1:numel(samples))
set(gca,'Yticklabel',{SampleNames(samples).Sample})
colormap([tamismagenta;tamisgreen]);



%compare two times points from same patient
samepatient=[2 3];
othersamples=[6 5 7];
samples2=[2 3 6 5 7];
figure(488); clf; hold on;
h1=bar([1 2], avgdistancefromPTanc(samepatient)*GenomeLength./sum(callablePos(:,samepatient)), .5);
h2=bar([3 4 5],avgdistancefromPTanc(othersamples)*GenomeLength./sum(callablePos(:,othersamples)), .5);

%xlabel('Average number of SNPs to patient LCA (?p)')
set(gca,'Xtick',1:numel(samples2))
set(gca,'Xticklabel',{'P2', 'P2T', 'P1','P3', 'P4'})
set(gca,'YAxisLocation','Right')
set(gca, 'Ytick', 0:5.25:21)
set(gca, 'Yticklabel', [0:5.25:21]/2.1)
xlabel('Average years to patient LCA (?p)')
set(h1,'FaceColor', rgb('Black'), 'EdgeColor', 'None')
set(h2,'FaceColor', rgb('Gray'), 'EdgeColor', 'None')
axis([0 numel(samples2)+1 0 25 ])
grid(gca)
set(gca, 'GridLineStyle', ':');




%% Display clickable scatter for data analysis of reverse strand versus forward strand

% for i=1:NSample
%     fwdMax=counts(sub2ind(size(counts),squeeze(maNT(:,i))',1:length(p),i*ones(length(p),1)'));
%     revMax=counts(sub2ind(size(counts),squeeze(maNT(:,i))'+4,1:length(p)',i*ones(length(p),1)'));
%     fwdCounts=sum(counts(1:4,:,i),1);
%     revCounts=sum(counts(5:8,:,i),1);
%     div_clickable_scatter_sigcolor(fwdMax./fwdCounts, revMax./revCounts, 'Forward strand MAF', 'Reverse strand MAF', i, parameters, counts, fwindows, positions, mutations, RefGenome, ScafNames, SampleNames, isolate_directory);
% end





%% Grid of mutations



%make matrix for easy comparison (better than above instances of the
%samefdn


figure(622); clf; hold on;

subm=zeros(numel(multiple_locations_within_sample_all),7);
polym=zeros(numel(multiple_locations_within_sample_all),7);
for i=1:numel(multiple_locations_within_sample_all)
    
    y=28-i;
    sub=sum(gnmat2==multiple_locations_within_sample_all(i) & fixedpos>0);
    poly=sum(gnmat2==multiple_locations_within_sample_all(i) & goodpos>0);
    
    f=find(gn==multiple_locations_within_sample_all(i),1);

    
    disp(genesizes(multiple_locations_within_sample_all(i)))
    for j=samples
        k=find(samples==j); 
        if sub(j)==0 & poly(j)==0 
            p1=patch([k; k+1; k+1; k], [y; y; y+1; y+1],rgb('White'));
            set(p1,'EdgeColor',rgb('Black'),'LineWidth',1)
        elseif sub(j)>0
            p1=patch([k; k+1; k+1; k], [y; y; y+1; y+1],tamisgreen);
            set(p1,'EdgeColor',rgb('Black'),'LineWidth',1)
        elseif poly(j)==1
            p1=patch([k; k+1; k+1; k], [y; y; y+1; y+1],tamismagenta);
            set(p1,'EdgeColor',rgb('Black'),'LineWidth',1)
        else       
            p1=patch([k; k+1; k+1; k], [y; y; y+1; y+1],rgb('DarkMagenta'));
            set(p1,'EdgeColor',rgb('Black'),'LineWidth',1)
            text(k+.3, y+.5, num2str(poly(j)), 'Color', rgb('White'), 'FontWeight', 'bold')
            disp(poly(j))
            %PLOT NUMBERS
          %  text()
        end
    end
    text(numel(samples)+1.2, y+.5, combined_annotation_all(f).protein )
    subm(i,:)=sub;
    polym(i,:)=poly;
    
end
axis([0 40 0 30])

[~,firstp]=ismember(multiple_locations_within_sample_all,gn);

multiple_annotations=combined_annotation_all(firstp);

save('formatrix', 'subm','polym', 'multiple_annotations')


% display information for annotating genes

for i=1:numel(multiple_locations_within_sample_all)
    j=find(gn==multiple_locations_within_sample_all(i),1);
     disp(combined_annotation_all(j).gene)
     disp(combined_annotation_all(j).protein)
%     disp([combined_annotation_all(gn==multiple_locations_within_sample_all(i)&sum(goodpos,2)>0).muts])
%     disp(genesizes(multiple_locations_within_sample_all(i)))
    disp(genesizes(multiple_locations_within_sample_all(i))/max(polym(i,:)))
end

% display fasta file

for i=1:numel(multiple_locations_within_sample_all)
    j=find(gn==multiple_locations_within_sample_all(i),1);
    disp(['>' combined_annotation_all(j).gene ]) %'-' combined_annotation_all(j).protein])
    disp(combined_annotation_all(j).translation)
end


%% Histograms of filtered data, histogram of isolates



goodposv=reshape(goodpos,numel(goodpos),1);
minorAFv=reshape(minorAF,numel(maf),1);
goodaf=minorAFv(goodposv>0);




%Plot
figure(21); clf;
k=1;
for i=samples(1:end-1);
    
    subplot(numel(samples),1,k); hold on;

   
    [freq, bins]=hist(expectedNumberMultipleMutatedGenes(:,i),max(expectedNumberMultipleMutatedGenes(:,i)));
    bar(bins-.5,freq,'FaceColor', rgb('gray'), 'EdgeColor', 'w')
    
    bar(numel(multiple_locations_within_sampleperlength{i}),numtrials,.25,'FaceColor',rgb('Red'),'EdgeColor', rgb('Red'))
    
    %title(SampleNames(i).Sample)
      axis([-0.6 8 -.01 numtrials])
     set(gca,'XTick', 0:2:8)
    set(gca,'YTick', [.25*numtrials .5*numtrials .75*numtrials numtrials])
    set(gca,'YTicklabel', {'.25', '.50', '.75', '1'})
    k=k+1;
end;


%all samples

figure(321); clf;
hist(goodaf,50);
xlabel('Minor allele frequency')
ylabel('Number of positions passing strict quality filters')
title('All samples')
axis([0 .5 -inf inf])


%each sample
for j=1:numel(samples)
    
    i=samples(j);
    disp(i)
    
    subplot(numel(samples),1,j); hold on;
    %figure(801);clf; hold on;
    %set(801,'Position',[100 scrsz(4)/2 scrsz(3)/5 scrsz(4)/6]);clf;hold on;

    b= hist([mutAF(goodpos(:,i)>0,i); 1.2],40);
    hist([mutAF(goodpos(:,i)>0,i); 1.2],40);
    h=findobj(gca,'Type','patch');
    set(h,'FaceColor', [1 0 1], 'EdgeColor', 'w')
    bar(.985,sum(fixedpos(:,i)),.03,'FaceColor',tamisgreen,'EdgeColor', 'w');
   % xlabel('Mutation allele frequency')
    %ylabel('Number of positions passing strict quality filters')
    %title(SampleNames(i).Sample)
    axis([0 1 -inf max(b)+5])
    disp(SampleNames(i).Sample)
    disp([sum(fixedpos(:,i)) sum(goodpos(:,i))])
    
end
set(gca,'Ytick',[15 30 45 60])

%isolates are still in terms of major allele frequency 
%isolates from sample 2
figure(850);clf;
hist(iMAF,10) 
xlabel('Major allele frequency')
ylabel('Number of positions')
title('Sputum2 - isolate estimate')
axis([.5 1 -inf inf])



%% Isolates divergence

isolateRef=repmat(NTs(ancnt(ismember(p,ip)))',1,numel(isolateNames));
isolatedistance=sum(isolateCalls~=isolateRef & isolateCalls~='N') *GenomeLength./ isolatesCoverage';

figure(870); clf; hold on;
hist([-8 isolatedistance],12)
set(gca,'Xtick',0:5:max(isolatedistance))
axis([0 max(isolatedistance)+1 -inf inf])
colormap([rgb('Gray')])
axis([0 max(isolatedistance)+1 -inf 14.5])
h=bar(mean(isolatedistance), 50, .5);
set(h, 'FaceColor', 'k', 'EdgeColor', 'None');
grid(gca)
set(gca, 'GridLineStyle', ':');
set(gca,'XAxisLocation','Top')
set(gca,'Xtick',0:2.1*2.5:31.5)
set(gca,'Xticklabel',0:2.5:15)


%% Simulation for overlap with previous study

figure(43); clf; hold on;
numtrials=1000;
overlap=zeros(numtrials,1);
%Pick 16 genes from genes mutated at least once
g=unique(gn(sum(goodpos(:,nonmutatorsamples),2)+sum(fixedpos(:,nonmutatorsamples),2)>0));
for i=1:numtrials
    x=g(floor(numel(g)*rand(numel(multiple_locations_within_sample_all),1)+1));
    overlap(i)=numel(intersect(prevmutated, x));
end
hist(overlap)


%repeat, removing patient 14/15/J from consideration
numtrials=1000;
overlap=zeros(numtrials,1);
%Pick 17 genes from genes mutated at least once
g=unique(gn(sum(goodpos(:,[6 5 7]),2)+sum(fixedpos(:,[6 5 7]),2)>0));
for i=1:numtrials
    x=g(floor(numel(g)*rand(numel(multiple_locations_within_sample_all)-2,1)+1)); %2 genes are specific to this patient
    overlap(i)=numel(intersect(prevmutated, x));
end




%%Simulation for overlap with multiple

figure(43); clf; hold on;
numtrials=1000;
overlap=zeros(numtrials,1);
%Pick numel(multiple) genes from genes mutated at least once
g=unique(gn(sum(goodpos(:,nonmutatorsamples),2)+sum(fixedpos(:,nonmutatorsamples),2)>0));
for i=1:numtrials
    x=g(floor(numel(g)*rand(numel(multiple),1)+1));
    overlap(i)=numel(intersect(multiple_locations_within_sample_all, x));
end
hist(overlap)
bar(numel(intersect(multiple, multiple_locations_within_sample_all)),numtrials, .25,'r');
    


%% Simulation for overlap with between patients


figure(44); clf; hold on;
numtrials=1000;
overlap2=zeros(numtrials,1);
%Pick numel(multiple) genes from genes mutated at least once
for i=1:numtrials
    x=g(floor(numel(g)*rand(numel(multiple),1)+1));
    overlap2(i)=numel(intersect(multiple_locations_within_sample_all, x));
end
hist(overlap2)
bar(numel(intersect(multiple, multiple_locations_within_sample_all)),numtrials, .25,'r');
    





%% Make table for adjaceny matrix and network

%cytoscape format

totalm=polym+subm;


for i=1:numel(multiple_locations_within_sample_all)
    k=1;
    for j=samples
        
        if subm(i,j)>0
            disp([num2str(multiple_locations_within_sample_all(i)) ' sub P' num2str(k)])
        elseif polym(i,j)==1
            disp([num2str(multiple_locations_within_sample_all(i)) ' poly P' num2str(k)])
        elseif polym(i,j)>1
            disp([num2str(multiple_locations_within_sample_all(i)) ' mpoly P' num2str(k)])
        end
         k=k+1;
    end
   
end




