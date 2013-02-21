function gn = bdolosa_genenums(positions,codingsequences, ChrStarts)

p=p2chrpos(positions,ChrStarts);


gn= zeros(size(p,1),1) ;
for c=1:numel(codingsequences)
    z = p(:,1)==c ;
    
    cgenes=genomic_position(codingsequences{c} ,p(z,2));
    
    %remove intragenic regions
    intragenic=cgenes~=floor(cgenes);
    x=find(z);  z(x(intragenic))=0;
    cgenes(intragenic)=[];
    
    cgenenums=char({codingsequences{c}(cgenes).gene}');
    gn(z) = str2num(cgenenums(:,6:end));
    
end
gn(gn==0)=[]; %remove empty regions

