function [df,MutDisMin,MutDisMax] = ana_strain_distance(Call,ttls)
NStrain = size(Call,2) ;
df = zeros(NStrain,NStrain) ;
for i=1:NStrain
    for j=1:NStrain
        df(i,j) = sum( Call(:,i)~=Call(:,j) & Call(:,i)~='N'  & Call(:,j)~='N' ) ;
    end
end

if nargout>1
    MutDisMin = zeros(size(Call,1),1) ;
    MutDisMax = zeros(size(Call,1),1) ;
    for k=1:size(Call,1)
        if (length(unique([Call(k,:), 'N']))<=2) ;
            MutDisMin(k) = nan ;
            MutDisMax(k) = nan ;
        else
            c=Call(k,:) ;
            c=c(ones(1,NStrain),:) ; d=c' ;
            MutDisMin(k) = min(df(c~=d & c~='N' & d~='N')) ;
            MutDisMax(k) = max(df(c~=d & c~='N' & d~='N')) ;
        end
    end
end


if nargin>1
    imagesc(log10(df))
    set(gca,'xtick',1:NStrain,'xticklabel',ttls) ;
    axis square
    colorbar
    caxis([0 3])
end

return