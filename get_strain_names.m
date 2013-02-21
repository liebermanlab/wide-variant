function StrainNames = get_strain_names

[~,~,tmp] = xlsread('Strain_Names.xls') ;

z = cellfun(@(V) any( isnan(V(:))), tmp(1:end,9) ) ;
z=~z ; z(1)=false ;

StrainNames.Roy_num = [tmp{z,1}] ;
StrainNames.Case = [tmp{z,2}] ;
StrainNames.StrainName = tmp(z,3) ;
StrainNames.PFGE = [tmp{z,4}] ;
StrainNames.SNPFileName = tmp(z,9) ;

return