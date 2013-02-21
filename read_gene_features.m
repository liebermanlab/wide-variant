function df = read_gene_features(f) 

for i=1:length(f)
    df(i,1) = break_str(f(i).Header) ;
    df(i,1).Sequence = f(i).Sequence ;
end

return

function features = break_str(str)

ftrs={
    '[gene='
    '] [protein='
    '] [protein_id='
    '] [location='
    ']'} ;
 
s=str ; 
x = 1 ; 
for i=1:length(ftrs)
    zz = min(findstr(str(x:end), ftrs{i})) ;
    x=x+zz+length(ftrs{i})-1 ;
    xx(i) = x ;
end
brk_s={} ;
for i=1:length(ftrs)-1
    brk_s{i} = str(xx(i):xx(i+1)-length(ftrs{i+1})-1) ;
end

features.gene = brk_s{1} ;
features.protein = brk_s{2} ;
features.protein_id = brk_s{3} ;
if brk_s{4}(1)=='c'       % if 'true', gene is found on complement strand. ('features.strand = false'). 
    loc = brk_s{4}(12:end-1) ;     % used to skip 'complement('
    features.strand=false ;         
else
    loc = brk_s{4} ;
    features.strand=true ;
end

features.loc1 = str2num(loc(1:find(loc=='.',1)-1)) ;
features.loc2 = str2num(loc(find(loc=='.',1)+2:end)) ;
features.Sequence =[] ;

return





