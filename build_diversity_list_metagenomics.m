

%% Set filtering parameters

loose_parameters=struct('minorfreqthreshold',.01, 'maxreads_perstrand',2000,...
    'minreads_perstrand',30,'min_bq',20,'min_mq' ,31, 'min_td', 8,...
    'max_td',42, 'max_sbp', 10,'max_bqp', 20,'max_mqp', 120,'max_tdp',15,...
    'max_percent_indels', .50);

strict_parameters=struct('minorfreqthreshold',.03, 'maxreads_perstrand',1300,...
    'minreads_perstrand',30,'min_bq',20,'min_mq' ,34, 'min_td', 8,...
    'max_td',42, 'max_sbp', 4,'max_bqp', 15,'max_mqp', 120,'max_tdp',15,...
    'max_percent_indels', .20);



%% Important variables

%Control is first sample in case


path('../scripts',path)
filen='diversity.mat';
RefGenome = 'Markers1';
processed=1;



ChrStarts=[0,3413074,5589203];
GenomeLength=6420401;



%% Get SampleNames and ScafNames 

SampleNames = read_sample_names ;
NSample = length(SampleNames) ;
RefGenome = {} ;
SampleDirs = {} ;

for i=1:NSample
    SampleDirs{i} = ['../' SampleNames(i).ExperimentFolder '/' SampleNames(i).Sample '/' SampleNames(i).AlignmentFolder ] ;
end

fr = fastaread(['../Reference_Genomes/' RefGenome '/genome.fasta']) ;

ScafNames = {fr.Header} ; 
for i=1:length(ScafNames) 
    f=find(ScafNames{i}==' ',1) ;
    ScafNames{i} = ScafNames{i}(1:f-1) ;
end




%% Build position list, display histrograms of each quality score

t=num2str(loose_parameters.MAFthreshold);(1:2)=[]; 
savename=['Positions_' t];
switch processed
    case 1
        [dp, positions, numfields] = find_diverse_positions_metagenomics(filen, loose_parameters, SampleDirs, SampleNames, ChrStarts);
        save(savename, 'p', 'positions', 'numfields')
    case 2
        load(savename)
end

fprintf(1,'Found %g diverse positions that meet specified criteria \n',length(dp)) ;


stop
%ADRIAN STOP HERE



%% Build counts list, annotate mutations

tic;
savename=['Counts_' t];

switch processed
    case 1
        [geneloc, cds, mutations] = annotate_mutations_auto_gb_metagenomics(positions,ScafNames,RefGenome) ;
        [counts, windows] = generate_diversity_struct_metagenomics(filen, SampleDirs, p, numfields, ChrStarts) ;
        save(savename,'counts', 'windows', 'mutations', 'geneloc', 'cds')
    case 2
        load(savename)
end
toc;



%% Calculate important things

[maf, maNT, minorNT] = div_major_allele_freq(counts);
goodpos=div_test_thresholds(counts,strict_parameters);

goodposv=reshape(goodpos,numel(goodpos),1);
mafv=reshape(maf,numel(maf),1);
goodmaf=mafv(goodposv>0);




%% Table
[combined_annotation_S, ~, ~, ~] =div_clickable_table_only_S(mutations, dp, ip, isolateCalls, counts,  windows, isolate_windows, strict_parameters, RefGenome, ScafNames, SampleNames, ChrStarts);

[combined_annotation_goodpos, combined_annotation_all, updated_annotation_diverse,  updated_annotation_isolates] =div_clickable_table(mutations, dp, ip, isolateCalls, counts,  windows, isolate_windows, strict_parameters, RefGenome, ScafNames, SampleNames, ChrStarts);



%% Create expected dNdS

%expected =div_calculate_expected(100, combined_annotation_goodpos, RefGenome, ScafNames, GenomeLength, ChrStarts);

expected=3.91;



%% Histograms of MAF as a function of N, S, or I


types=[combined_annotation_all.type];
itypes=[updated_annotation_isolates.type];


step=.025;
edges=.5:step:1;
allN=zeros(numel(edges),1);
allS=zeros(numel(edges),1);
allI=zeros(numel(edges),1);
allP=zeros(numel(edges),1);

%each sample & all
for i=2:NSample
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
    legend('N', 'S','I','P', 'Location', 'NorthWest')
    
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
title('All Samples')
axis([.5 1 -inf inf])
legend('N', 'S','I','P', 'Location', 'NorthWest')


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
legend('N', 'S','I','P', 'Location', 'NorthWest')


%% dN/dS as a function of MAF cutoff


%there is a slight error as polymorphism could be S in one sample and N in
%another-- recorded as S here.

Nbins=20;
dNdS=nan(Nbins,NSample);
N=nan(Nbins,NSample);
S=nan(Nbins,NSample);

f=.5+(1:Nbins)*(.5/Nbins);
for i=2:NSample
    for j=1:Nbins
        N(j,i)=sum(types(maf(:,i)<f(j) & goodpos(:,i) > 0)=='N');
        S(j,i)=sum(types(maf(:,i)<f(j) & goodpos(:,i) > 0)=='S');
        if S(j,i) > 0 && N(j,i) > 0
            dNdS(j,i)=(N(j,i)/S(j,i))/expected;
        end
    end
end
figure(400);clf;hold on;
plot(f,dNdS)
legend({SampleNames.Sample})
xlabel('Major allele frequency cutoff (positions with MAF <x considered)')
ylabel('dN/dS (normalized to expectation)')


%expected = .7963

%all = .71-.79
%<.75 = .70 -.86
%all, w/o sample 3 = .76-.85
%<.75, w/o sample 3 = .76-.94



%over all samples
N(isnan(N))=0;
S(isnan(S))=0;
allN=sum(N,2);
allS=sum(S,2);
alldNdS=nan(Nbins,1);
for i=1:numel(allN)
    if allS(i) > 0 
        alldNdS(i)=(allN(i)./allS(i))/expected;
    end
end
figure(860);plot(f,alldNdS);


%dNdS from isolates
idNdS=nan(Nbins,1);

for j=1:Nbins
    N=sum(itypes(iMAF<f(j))=='N');
    S=sum(itypes(iMAF<f(j))=='S');
    if S > 0 && N > 0
        idNdS(j)=(N/S)/expected;
    end
end
figure(410);clf;hold on;
plot(f,idNdS)
legend('isolates')
xlabel('Major allele frequency cutoff (positions with MAF <x considered)')
ylabel('dN/dS (normalized to expectation)')




%% Are N and S in different places in genes?

figure(723); clf; hold on;
Ndist=[];
Sdist=[];
for i=2:NSample
    Ndist=[Ndist [combined_annotation_all(types(maf(:,i)<f(j) & goodpos(:,i) > 0)=='N').nt_pos]];
    Sdist=[Sdist [combined_annotation_all(types(maf(:,i)<f(j) & goodpos(:,i) > 0)=='S').nt_pos]];
end

hist(Ndist,50)
h=findobj(gca,'Type','patch');
set(h,'FaceColor', 'r', 'EdgeColor', 'w')
hist(Sdist,50)
legend({'N' 'S'})
xlabel('nt pos in gene')
ylabel('times found')



%% Display clickable scatter for data analysis of reverse strand versus forward strand

for i=1:NSample
    fwdMax=counts(sub2ind(size(counts),squeeze(maNT(:,i))',1:length(p),i*ones(length(p),1)'));
    revMax=counts(sub2ind(size(counts),squeeze(maNT(:,i))'+4,1:length(p)',i*ones(length(p),1)'));
    fwdCounts=sum(counts(1:4,:,i),1);
    revCounts=sum(counts(5:8,:,i),1);
    div_clickable_scatter_sigcolor__metagenomics(fwdMax./fwdCounts, revMax./revCounts, 'Forward strand MAF', 'Reverse strand MAF', i, strict_parameters, counts, windows, positions, mutations, RefGenome, ScafNames, ChrStarts, SampleNames);   
end









