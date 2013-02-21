function expected = div_calculate_expected(Ntimes, actual_mutations, Ref, Scaf, GLength, ChrStarts)

% May 2012, Tami Lieberman


Npositions=numel(actual_mutations)*Ntimes;
promoterdistance=150;

%generate random positions on genomes
p=floor(rand(Npositions,1)*GLength)+1;
positions=p2chrpos(p,ChrStarts);


%use actual nts, use each pair Ntimes, must restructure
nts={actual_mutations.nts};
n=cellfun(@numel,nts); for i=1:find(n>2); nts{i}=nts{i}(1:2); end;    
nts=double(char(nts));
nts(nts==65)=1; nts(nts==84)=2; nts(nts==67)=3; nts(nts==71)=4;

nt_i=mod(1:Npositions,numel(actual_mutations))+1;

%finds which genes these random mutations are in
[~, ~, annotations] = annotate_mutations_auto_gb(positions,Scaf,Ref);



%count types of mutations
N=0;
S=0;
P=0;
I=0;
U=0;



for i=1:Npositions
        
    %intragenic
    if floor(annotations(i).gene_num)==annotations(i).gene_num
        if numel(annotations(i).AA)==4
            annotations(i).annotation=annotations(i).protein;
            if strcmp(annotations(i).AA(nts(nt_i,1)),annotations(i).AA(nts(nt_i,2)));
                S=S+1;
            else
                N=N+1;
            end
        else
            U=U+1;
        end
        
    %intergenic
    else
        if (~isempty(annotations(i).distance1) && annotations(i).distance1 > -1*promoterdistance && annotations(i).distance1 < 0) || ...
                (~isempty(annotations(i).distance2) && annotations(i).distance2 > -1*promoterdistance && annotations(i).distance2 < 0)
            P=P+1;
        else
            I=I+1;
        end
    end
    
end


expected=N/S;


end
