function par = read_vcf_line(lin)

par = VCF_struct ;

f = [0, find(lin==9), length(lin)+1];

str={} ;
for i=1:length(f)-1
    str{i} = lin(f(i)+1:f(i+1)-1) ;
    if strcmp(str{i},'.'), str{i}=nan ; end
end

par(1).scaf = str{1} ;
par.pos = str2num(str{2}) ;
par.id = str{3} ;
par.ref = str{4} ;
par.alt = str{5} ;
par.qual = str2num(str{6}) ;
par.filter = str{7} ;
par.info = str{8} ;
t = read_vcf_info(str{8}) ;
f = fieldnames(t) ;
for i=1:length(f)
    par.(f{i}) = t.(f{i}) ;
end


%TDL 2012 unclear to me what this gen1 and gen2 is supposed to be it seems
%to be a 0/1 of whether or not there was information at this position?
%gen2 hold the info in 'GT'
if length(str)>=10
    par.gen1 = str{9} ;
    if strcmp(par.gen1,'PL')
        %***
        par.gen2 = 1 ; 
    else
        par.gen2 = str2num(str{10}(1:find(str{10}==':',1)-1)) ;
    end
end
par.lins = [] ; % lin ;
end