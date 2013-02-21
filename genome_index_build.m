function genome_index_build(fldr)

if nargin<1
   fldr = '.' ;
end

cdr = pwd ;

cd(fldr)

if exist('genome.fasta','file')
    delete('genome.fasta')
end

fastanames = dir('*.fasta') ;

fastanames = { fastanames.name } ;

cat_line = '' ; 
for i=1:length(fastanames)
    cat_line = [cat_line ' '  fastanames{i}] ;  
end
cat_line(1) = [] ; 

eval (['!cat ' cat_line ' >  genome.fasta'])
    
!/opt/bowtie/bowtie-build genome.fasta genome
!/opt/bowtie2/bowtie2-build genome.fasta genome


cd(cdr) ;

end