function filter_reads(fname_in,fname_out,Phred_offset,method,param)



disp('runningfilterreads')
tic ;

fin = fopen(fname_in,'r') ;
fout = fopen(fname_out,'w') ;
flog = fopen([fname_out '.log'],'w') ;

y = 0 ;
w = 0 ;

while ~feof (fin)
    w=w+1 ;
    if ~mod(w,50000), fprintf(flog,'.'); end
    t = fgets(fin) ;
    s = fgets(fin) ; %actual sequence
    p = fgets(fin) ;
    q = fgets(fin) ; %phred scores
    switch method
        case 'all_quality' % Phred = param(1)
            if all(q(1:end-1)-Phred_offset >= param(1))
                y = y + 1 ;
                fwrite(fout,[t s p q]) ;
            end
        case 'trim_end' % Phred=param(1), Min Length = param(2)
            f = find(q(1:end-1) < param(1)+ Phred_offset,1) ;
            if isempty(f), f=length(q) ; end
            if f>param(2)
                y = y + 1 ;
                fwrite(fout,[t s([1:f-1,end]) p q([1:f-1,end])]) ;
            end
    end
end

fclose(fin) ;
fclose(fout) ;

fprintf(flog,'\n');
fprintf(flog, ['Method ' method '\nParam: ' num2str(param) '\nNumber of lines: %g / %g; Time: %5.0f min \n'], y, w, toc/60) ;

fclose(flog) ;
end


