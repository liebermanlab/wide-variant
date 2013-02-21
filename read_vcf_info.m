function info = read_vcf_info(vcf_info) ;

%DP=270;AF1=0;AC1=0;DP4=167,101,0,2;MQ=20;FQ=-282;PV4=0.14,2e-14,1,0.21

vcf_info = [';' vcf_info ';' ];
f1 = find(vcf_info=='=') ;
f2 = find(vcf_info==';') ;

info.DP=nan ;
info.AF1=nan ;
info.AC1=nan ;
info.DP4=[nan,nan,nan,nan] ;
info.MQ=nan ;
info.FQ=nan ;
info.PV4=[nan,nan,nan,nan] ;

for i=1:length(f1)
    nm = vcf_info(f2(find(f2<f1(i),1,'last'))+1:f1(i)-1) ;
    info.(nm) = str2num(vcf_info(f1(i)+1:f2(find(f2>f1(i),1))-1)) ;
end

end