function AlignmentTable = read_alignment_table

d = tdfreadunix('alignment_params.csv',',') ;

%***
%if ~all(ismember({'Batch','Lane','Barcode','Sample','Alignments'},fieldnames(d)))
%    error('Notice headers in samples.csv')
%end

k = {} ;
for f=fieldnames(d)'
    if ischar(d.(f{1}))
        d.(f{1}) = cellstr(d.(f{1})) ;
    else
        d.(f{1}) = num2cell(d.(f{1})) ;
    end
    k(end+1,:) = {d.(f{1}){:}} ;
end
AlignmentTable = cell2struct(k,fieldnames(d),1) ;



end