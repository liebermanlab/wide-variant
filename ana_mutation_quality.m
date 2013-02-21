function MutQual = ana_mutation_quality(Call,Qual, Fig, strains)

%November 2012
%added variable strains so that an outgroup can be removed


if nargin<3
    Fig = 0 ;
end

NStrain = size(Call,2) ;

if nargin<4
    strains = 1:NStrain ;
end


MutQual = zeros(size(Call,1),1) ; 
for k=1:size(Call,1)
    if (length(unique([Call(k,strains), 'N']))<=2) ; %if there is only one type of non-N call, skip this location
        MutQual(k) = nan ;
    else
        c=Call(k,:) ; c1=c(ones(1,NStrain),:) ; c2=c1' ;
        q=Qual(k,:) ; q1=q(ones(1,NStrain),:) ; q2=q1' ;
        g=c1~=c2 & c1~='N' & c2~='N' ; %find places where calls disagree
        MutQual(k) = max(min(q1(g),q2(g))) ;%min(q1(g),g2(g)) gives lower qual for each disagreeing pair of calls, we then find the best of these
    end
end

MutQual(isnan(MutQual))=0;


if Fig
    figure(Fig); clf;
    semilogx(sort(MutQual),length(MutQual):-1:1,'.') ;
    hold on
    xlabel('Mutation Quality')
    ylabel('Number of Mutations')
    axis([0 max(MutQual)*1.05 0 1.05*sum(MutQual>0)])
    grid
end

return



