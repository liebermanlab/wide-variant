function lins = read_vcf_file(filename)
Tab=9 ;
y=0 ;

fid = fopen(filename,'r') ;



i=0 ;
while ~feof(fid) 
    s = fgetl(fid) ;
    f = find(s==Tab,2) ;
    if length(f)>=2
        if i>0
            lins(i) = read_vcf_line(s) ;
        end
        i=i+1 ;
    end
end

if ~exist('lins', 'var')
      lins=0;
end
    
fclose(fid) ;

end

