function par = read_vcf_line_light(lin)

par = VCF_struct_light ;

f = [0, find(lin==9), length(lin)+1];

str={} ;
for i=1:length(f)-1
    str{i} = lin(f(i)+1:f(i+1)-1) ;
    if strcmp(str{i},'.'), str{i}=nan ; end
end

par(1).scaf = str{1} ;
par.pos = str2num(str{2}) ;
par.alt = str{5} ;
par.ref = str{4} ;
t = read_vcf_info(str{8}) ;
par.FQ = t.FQ ;


end