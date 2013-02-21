%% Load mutgen
%***
%load mutgenvcf_3

%% quality for each mutation

MutQual = zeros(size(Call3,1),1) ; 
for k=1:size(Call3,1)
    if (length(unique([Call3(k,:), 'N']))<=2) ;
        MutQual(k) = nan ;
    else
        c=Call3(k,:) ; c1=c(ones(1,NStrain),:) ; c2=c1' ;
        q=Qual3(k,:) ; q1=q(ones(1,NStrain),:) ; q2=q1' ;
        g=c1~=c2 & c1~='N' & c2~='N' ;
        MutQual(k) = max(min(q1(g),q2(g))) ;
    end
end

figure(6);clf
semilogx(sort(MutQual),length(MutQual):-1:1,'.') ;
hold on
xlabel('Mutation Quality')
ylabel('Number of Mutations')
axis tight
grid

%% quality threshold analysis (movie)
qual_th = 0:20:200 ; 
df_m = zeros(NStrain,NStrain,length(qual_th)) ;
for i=1:NStrain
    for j=1:NStrain
        for m=1:length(qual_th)
            df_m(i,j,m) = sum( Call3(:,i)~=Call3(:,j) & Call3(:,i)~='N'  & Call3(:,j)~='N' & Qual3(:,i)>qual_th(m) & Qual3(:,j)>qual_th(m)) ;
        end
    end
end
figure(7);clf
set(gcf,'position',[10   570   560   420])
z=1:Nstrain ;
for m=1:size(df_m,3)
    imagesc(log10(df_m(z,z,m)))
    title(sprintf('Quality > %g',qual_th(m))) ;    
    set(gca,'xtick',1:NStrain,'xticklabel',{StrainNames(z).Sample}) ;
    axis square
    colorbar
    caxis([0 3])
    %pause
end

%% mutation list

z=[1 2 3 4 7 11 12 15 17 18 19 20 21 22 9 8 16 13 14 5 6 10] ;  %***
qual_0 = 35 ; %***

Call4=Call3 ;
Call4(Qual3<qual_0) = 'N' ; 

df = zeros(NStrain,NStrain) ;
for i=1:NStrain
    for j=1:NStrain
        df(i,j) = sum( Call4(:,i)~=Call4(:,j) & Call4(:,i)~='N'  & Call4(:,j)~='N' ) ;
    end
end

figure(8);clf
set(gcf,'position',[10   570   560   420])
imagesc(df(z,z))
title(sprintf('Quality > %g',qual_0)) ;
set(gca,'xtick',1:NStrain,'xticklabel',{StrainNames(z).Sample}) ;
axis square
colorbar

MutDisMin = zeros(size(Call3,1),1) ; 
MutDisMax = zeros(size(Call3,1),1) ; 
for k=1:size(Call4,1)
    if (length(unique([Call4(k,:), 'N']))<=2) ;
        MutDisMin(k) = nan ;
        MutDisMax(k) = nan ;
    else
        c=Call4(k,:) ;
        c=c(ones(1,NStrain),:) ; d=c' ;
        MutDisMin(k) = min(df(c~=d & c~='N' & d~='N')) ;
        MutDisMax(k) = max(df(c~=d & c~='N' & d~='N')) ;
    end
end

figure(9);clf
set(gcf,'position',[10   60   560   420])
hist(MutDisMin,16) ;
hold on
%***
MutDis_Th = 50 ;
plot(MutDis_Th*[1 1],ylim,'--k')
xlabel('Minimal mutation difference between strains')
ylabel('Number of mutations')
%gn =  [Mut3(~isnan(MutDisMin)).gene_num] ;
%title(sprintf('Qual>%g, Number of mutations: %g, Number of genes: %g',qual_0, sum(~isnan(MutDisMin)), length(unique(gn(gn==round(gn))))))
figure(6)
plot(qual_0,sum(~isnan(MutDisMin)),'or')

[~,ind]=sort(MutDisMin) ;
ind = ind(MutDisMin(ind)<MutDis_Th);

mutgen = {} ;
for i=1:length(z)
    % mutgen(1,1:length(z)) = [StrainNames.SNPFileName(z)] ;
    mutgen{1,i} = StrainNames(z(i)).Sample(end-2:end) ;
end
mutgen(2:length(ind)+1,1:length(z)) = num2cell(Call3(ind,z)) ;
mutgen{1,length(z)+1} = 'Gene' ;
mutgen(2:length(ind)+1,length(z)+1) = {Mut3(ind).gene} ;         
mutgen{1,length(z)+2} = 'Protein' ;
mutgen(2:length(ind)+1,length(z)+2) = {Mut3(ind).protein} ;
mutgen{1,length(z)+3} = 'NonSyn' ;
mutgen(2:length(ind)+1,length(z)+3) = {Mut3(ind).NonSyn} ;
mutgen{1,length(z)+4} = 'MutDisMin' ;
mutgen(2:length(ind)+1,length(z)+4) = num2cell(MutDisMin(ind)) ;
mutgen{1,length(z)+5} = 'MutQual' ;
mutgen(2:length(ind)+1,length(z)+5) = num2cell(MutQual(ind)) ;
mutgen{1,length(z)+6} = 'Scaf' ;
mutgen(2:length(ind)+1,length(z)+6) = {Mut3(ind).scafold} ;
mutgen{1,length(z)+7} = 'Pos' ;
mutgen(2:length(ind)+1,length(z)+7) = {Mut3(ind).pos} ;

% open mutgen
% xlswrite('Mutation_list',mutgen)

%%
figure(11);
MutGenVCF3ind = MutGenVCF3 ; 
DP4 = reshape([MutGenVCF3ind.DP4],4,size(MutGenVCF3ind,1),size(MutGenVCF3ind,2)) ;

set(gcf, 'Position',[10         50        1450         850]);
columnname =   mutgen(1,1:end);

%columnformat = {'numeric', 'bank', 'logical', {'Fixed' 'Adjustable'}};
%columneditable =  [false false true true]; 
t = uitable('Units','normalized','Position',[0.05 0.05 0.9 0.9], 'Data', mutgen(2:end,:),... 
            'ColumnName', columnname,...
            'RowName',[],'ColumnWidth',num2cell([ones(1,length(z))*30, 100, 200, 20, 20 ,20, 20, 20]), ...
        'CellSelectionCallback',@igv_caller );
%            'ColumnFormat', columnformat,...
%            'ColumnEditable', columneditable,...
set(t,'userdata',{[StrainNames.SNPFileName(z)],RefGenome,ScafNames,0,DP4(:,:,z)})


MutGenVCF3ind = MutGenVCF3(ind,:) 

for i=1:length(MutGenVCF3ind(:))
     t = read_vcf_info(MutGenVCF3ind(i).info) ;
    f = fieldnames(t) ;
    for j=1:length(f)
        MutGenVCF3ind(i).(f{j}) = t.(f{j}) ;
    end
end






