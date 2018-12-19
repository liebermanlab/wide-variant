function lins = read_vcf_file_dindel(filename)
Tab=9 ;
y=0 ;

fid = fopen(filename,'r') ;



i=0 ;
while ~feof(fid) 
    s = fgetl(fid) ;
    f = find(s==Tab,2) ;
    if length(f)>=2
        %if indel length is too long, it reads onto next line instead of
        %next tab
        if length(s) < 40
            s= [s fgetl(fid)] ;
        end
        if i>0
            lins(i) = read_vcf_line_dindel(s) ;
        end
        i=i+1 ;
    end
end

if ~exist('lins', 'var')
    lins=0;
end
    
fclose(fid) ;

end

