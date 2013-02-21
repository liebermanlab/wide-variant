function genes=locustags2numbers(locustags)

genes=zeros(numel(locustags),1);
for i=1:numel(genes);
    gn=locustags{i};
    if numel(gn)>6
        genes(i)=str2double(gn(end-4:end));
    else
        genes(i)=0;
    end
end
