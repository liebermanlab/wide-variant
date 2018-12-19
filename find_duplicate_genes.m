function [numbertimesfound, uniques] = find_duplicates(initiallist)

%works for operons or genes

uniques=unique(initiallist);
numbertimesfound=zeros(1,numel(uniques));
for k=1:numel(uniques)
    if isnumeric(uniques{1})
        numbertimesfound(k)=sum(uniques==initiallist(k));
    else
        numbertimesfound(k)=sum(strcmp(uniques{k},initiallist));
    end
end