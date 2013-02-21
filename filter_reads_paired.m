function filter_reads_paired(fname_in1,fname_in2,fname_out1,fname_out2,Phred_offset,method,param)



disp('runningfilterreads')
tic ;

fin1 = fopen(fname_in1,'r') ;
fin2 = fopen(fname_in2,'r') ;
fout1 = fopen(fname_out1,'w') ;
fout2 = fopen(fname_out2,'w') ;
flog = fopen([fname_out1 '.log'],'w') ;

y = 0 ;
w = 0 ;

while ~feof (fin1)
    w=w+1 ;
    if ~mod(w,50000), fprintf(flog,'.'); end
    t1 = fgets(fin1) ;
    s1 = fgets(fin1) ; %actual sequence
    p1 = fgets(fin1) ;
    q1 = fgets(fin1) ; %phred scores
    t2 = fgets(fin2) ;
    s2 = fgets(fin2) ; %actual sequence
    p2 = fgets(fin2) ;
    q2 = fgets(fin2) ; %phred scores
    switch method
        case 'all_quality' % Phred = param(1)
            if all(q1(1:end-1)-Phred_offset >= param(1)) && all(q2(1:end-1)-Phred_offset >= param(1))
                y = y + 1 ;
                fwrite(fout1,[t1 s1 p1 q1]) ;
                fwrite(fout2,[t2 s2 p2 q2]) ;
            end
        case 'trim_end' % Phred=param(1), Min Length = param(2)
            f1 = find(q1(1:end-1) < param(1)+ Phred_offset,1) ;
            f2 = find(q2(1:end-1) < param(1)+ Phred_offset,1) ;
            if isempty(f1), f1=length(q1) ; end
            if isempty(f2), f2=length(q2) ; end
            if f1>param(2) && f2>param(2)
                y = y + 1 ;
                fwrite(fout1,[t1 s1([1:f1-1,end]) p1 q1([1:f1-1,end])]) ;
                fwrite(fout2,[t2 s2([1:f2-1,end]) p2 q2([1:f2-1,end])]) ;
            end
    end
end

fclose(fin1) ;
fclose(fin2) ;
fclose(fout1) ;
fclose(fout2) ;

fprintf(flog,'\n');
fprintf(flog, ['Method ' method '\nParam: ' num2str(param) '\nNumber of lines: %g / %g; Time: %5.0f min \n'], y, w, toc/60) ;

fclose(flog) ;
end