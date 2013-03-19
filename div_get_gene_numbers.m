function genesN = div_get_gene_numbers(genestructure)

%Tami Lieberman
%May 2012

genesN=zeros(length(genestructure),1);

tags={genestructure.locustag}';
num_starts=cellfun(@locustag2num,tags);
for i = 1:length(genestructure)
    if num_starts(i)>0
        genesN(i,1)=str2double(tags{i}(num_starts(i):end));
    end
end

%genesN=unique(genesN);  %not sure why this was here previously, commented
%out
%3/6/13

end