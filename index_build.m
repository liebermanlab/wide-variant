function index_build


[~,~,d] = xlsread ('NCBI_Genomes.xls', 'Sheet1', '', 'basic') ; 

bowtie_line = '' ; 
for i=1:length(d)
    bowtie_line = [bowtie_line ' /media/my_passport/PhD/Automated/NCBI_Genomes/'  d{i} '.fasta'] ;  
end
bowtie_line(1) = [] ; 

eval (['!cat ' bowtie_line ' >  genome.fasta'])
    
!bowtie-build genome.fasta genome

end