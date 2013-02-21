function FilterTable = read_filter_table

d = tdfreadunix('read_filter_params.csv',',') ;

if ~all(ismember({'Filter','Method'},fieldnames(d)))
    error('Notice headers in read_filter_params.csv')
end

d.Params = [d.Param1, d.Param2, d.Param3] ;

d = rmfield(d,{'Param1','Param2','Param3'}) ;

k = {} ;

d.Filter = cellstr(d.Filter) ;
d.Method = cellstr(d.Method) ;
d.Params = num2cell(d.Params,2) ;
for f=fieldnames(d)'
    k(end+1,:) = {d.(f{1}){:}} ;
end
FilterTable = cell2struct(k,fieldnames(d),1) ;

end

