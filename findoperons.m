function operonnumbers=findoperons(ps,genomicoperons)


operonnumbers=zeros(size(ps));
for i=1:numel(ps);
    p=ps(i);
    if(find(genomicoperons(:,1)<= p & genomicoperons(:,2)>=p));
        possibleoperons=find(genomicoperons(:,1)<= p & genomicoperons(:,2)>=p);
        
        operon=1;
        
        %if overlapping operons, pick larger one
        if numel(possibleoperons)>1
            operonlengths=genomicoperons(possibleoperons,2)-genomicoperons(possibleoperons,1);
            operonnumbers(i)=possibleoperons(find(operonlengths==max(operonlengths),1));  %largest operon
        else
            operonnumbers(i)=possibleoperons;
        end
    end
end