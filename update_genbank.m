gb = fopen('Reference_Genomes/SaureusNCTC8325/sequence.gb','a');
seq = fastaread('Reference_Genomes/SaureusNCTC8325/genome.fasta');
len = length(seq.Sequence);

for i=1:60:len-60
    linenum = sprintf('%9s',num2str(i));
    fprintf(gb,[linenum ' ' lower(seq.Sequence(i:i+9)) ...
        ' ' lower(seq.Sequence(i+10:i+19)) ' ' lower(seq.Sequence(i+20:i+29)) ...
        ' ' lower(seq.Sequence(i+30:i+39)) ' ' lower(seq.Sequence(i+40:i+49)) ...
        ' ' lower(seq.Sequence(i+50:i+59)) '\n']);
end

i = i+60;
linenum = sprintf('%9s',num2str(i));
fprintf(gb,[linenum ' ']);

for j = i:10:len-10
    fprintf(gb,[lower(seq.Sequence(j:j+9)) ' ']);
end


for k = j+10:len
    fprintf(gb,lower(seq.Sequence(k)));
end

fprintf(gb,'\n');
fclose(gb);