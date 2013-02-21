global cMG cM

%% Load diversity data
%***
load diversitytable

numSamples=size(counts,1)
z=1:numSamples ;

%% quality for each mutation
figure(6);clf
MutQual = ana_mutation_quality(Call3,Qual3,1) ;

%% quality threshold analysis (movie)
qual_th = 0:20:200 ;
figure(7);clf
set(gcf,'position',[10   570   560   420])
for m=1:length(qual_th)
    Call4 = Call3 ; Call4(Qual3<qual_th(m)) = 'N' ;
    df = ana_strain_distance(Call4(:,z),{StrainNames(z).Sample}) ;
    title(sprintf('Quality > %g',qual_th(m))) ;
    %    pause
end

%% mutation list

qual_0 = 35 ; %***

Call4=Call3 ; Call4(Qual3<qual_0) = 'N' ;
[df,MutDisMin,MutDisMax] = ana_strain_distance(Call4) ;

figure(8);clf
set(gcf,'position',[10   570   560   420])
ana_strain_distance(Call4(:,z),{StrainNames(z).Sample}) ;
title(sprintf('Quality > %g',qual_0)) ;

figure(9);clf
set(gcf,'position',[10   60   560   420])
hist(MutDisMin,16) ;
hold on
%***
MutDis_Th = 50 ;
plot(MutDis_Th*[1 1],ylim,'--k')
xlabel('Minimal mutation difference between strains')
ylabel('Number of mutations')
gn = [Mut3(~isnan(MutDisMin)).gene_num] ;
title(sprintf('Qual>%g, Number of mutations: %g, Number of genes: %g',...
    qual_0, sum(~isnan(MutDisMin)), length(unique(gn(gn==round(gn))))))


figure(6)
plot(qual_0,sum(~isnan(MutDisMin)),'or')

[~,ind]=sort(MutDisMin) ;
ind = ind(MutDisMin(ind)<MutDis_Th);

Z = z ;

ana_clickable_table(StrainNames, Call3, Mut3, MutGenVCF3, ScafNames, RefGenome, ind, z,...
    {{'MutDisMin',num2cell(MutDisMin)},{'MutQual',num2cell(MutQual)}},11)



