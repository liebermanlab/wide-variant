function [call, qual] = vcfline2callqual(lin)


f = [0, find(lin==9), length(lin)+1];

str={} ;
for i=1:length(f)-1
    str{i} = lin(f(i)+1:f(i+1)-1) ;
    if strcmp(str{i},'.')
        str{i}=nan;
    end
end

alt=str{5};
ref=str{4};

if isempty(alt) %if no alt
    call = 'N' ;
elseif isnan(alt)
    if numel(ref)==1 & isempty(findstr(str{8},'INDEL')) %not indels
        call = ref;
    else
        call = 'D';
    end
elseif any(alt==',') %ambigious-- give N
    call = 'N' ;
elseif length(alt)==length(ref)
    call = alt ;
elseif length(alt)>length(ref)
    call = 'I' ;
else
    call = 'D' ;
end

t = read_vcf_info(str{8}) ;
qual = t.FQ ;


end