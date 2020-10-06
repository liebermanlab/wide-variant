function mut_annotations = annotate_mutations_gb(Positions,REFGENOMEFOLDER)


%% Written and modified by Idan Yelin, Roy Kishony, Tami Lieberman, and Hattie Chung

% May need revisions for your reference genome, especially subfunction read_gb_modified



%% Intialize



nts='atcg';
rc='tagc';
[ChrStarts, GLength, ~, ScafNames]= genomestats(REFGENOMEFOLDER);


if ~exist([REFGENOMEFOLDER '/cds_sorted.mat'], 'file')
    read_gb_modified(REFGENOMEFOLDER)
end

load([REFGENOMEFOLDER '/cds_sorted.mat'])

%%  Fill in annotation for mutated positions

%find corresponding genes by using reference created in creation of
%cds_sorted
p=chrpos2index(Positions,ChrStarts);
genenum=genomic_position_all(CDS, GLength, ChrStarts);
mut_genenum=genenum(p);



%fill in details
mut_annotations=[] ;
for i=1:size(Positions,1)
    
    
    mut_annotations(i).gene_num = mut_genenum(i) ;
    mut_annotations(i).scaffold = Positions(i,1) ;
    mut_annotations(i).pos = Positions(i,2) ;
    
    nScf = Positions(i,1) ;
    if mut_genenum(i)==round(mut_genenum(i)) % intragenic
      %  disp(i)
        cdf = CDS{nScf}(mut_genenum(i)) ;
        
        %if product is more than one line, reshape
        if size(cdf.product,1)>1
            cdf.product=reshape(cdf.product',size(cdf.product,1)*size(cdf.product,2),1)';
        end
        
        if isfield(cdf,'gene')
            mut_annotations(i).gene       = cdf.gene ;
        else
            mut_annotations(i).gene = '';
        end
        if isfield(cdf,'product')
            mut_annotations(i).protein    = cdf.product ;
        else
            mut_annotations(i).protein = '';
        end
        mut_annotations(i).protein_id = cdf.protein_id ;
        mut_annotations(i).strand     = cdf.strand ; %1 indicates reverse strand, 0 forward strand
        mut_annotations(i).loc1       = cdf.loc1 ;
        mut_annotations(i).loc2       = cdf.loc2 ;
        mut_annotations(i).Sequence   = cdf.Sequence ;
        mut_annotations(i).note       = cdf.note ;
        mut_annotations(i).locustag   = cdf.locustag ;
        if isfield(cdf,'oldlocustag')
            mut_annotations(i).oldlocustag   = cdf.oldlocustag ;
        end
        mut_annotations(i).translation = nt2aa(cdf.Sequence, 'GENETICCODE', 11,'ACGTOnly','F');
        
        
        
        if mut_annotations(i).strand
            p = double(mut_annotations(i).loc2) - double(Positions(i,2)) + 1;
        else
            p = double(Positions(i,2)) -double(mut_annotations(i).loc1) + 1;
        end
        mut_annotations(i).nt_pos = p ;
        
        aan = floor((p-1)/3) + 1 ;
        ncn = p-(aan-1)*3 ;
        codons=cell(4,1);
        AA='';
        mut_annotations(i).aa_pos=aan;
        
        if numel(mut_annotations(i).Sequence) >= aan*3 & mut_annotations(i).translation > 1;
            codon = mut_annotations(i).Sequence(aan*3-2:aan*3) ;
            for j=1:4 %for each nts
                if mut_annotations(i).strand
                    codon(ncn)=rc(j);
                else
                    codon(ncn)=nts(j);
                end
                codons{j}=codon;
                AA(j) = nt2aa(codon, 'ACGTOnly', false, 'ALTERNATIVESTARTCODONS','F', 'GENETICCODE', 11) ;
                %setting ACGTonly to false prevents nt2aa from throwing an error for non-ACTG calls
                %setting ALTERNATIVESTARTCODONS to false prevents nts from
                %being called as a start codon all the time
            end
        end
        
        
        mut_annotations(i).codons   = codons ;
        mut_annotations(i).AA   = AA ;
        mut_annotations(i).NonSyn = length(unique(AA))>1 ;
    else %intergenic
        if floor(mut_genenum(i))>=1 % gene before position
            cdf = CDS{nScf}(floor(mut_genenum(i))) ;
            if isfield(cdf,'gene')
                mut_annotations(i).gene1       = cdf.gene ;
            else
                mut_annotations(i).gene1 = '';
            end
            mut_annotations(i).locustag1 = cdf.locustag ;
            mut_annotations(i).distance1 = Positions(i,2) - cdf.loc2;
            mut_annotations(i).protein1    = cdf.product ;
            if cdf.strand==1
                mut_annotations(i).distance1 = mut_annotations(i).distance1 * -1;
            end
        end
        if ceil(mut_genenum(i))<=length(CDS{Positions(i,1)}) % gene after position
            cdf = CDS{nScf}(ceil(mut_genenum(i))) ;
            if isfield(cdf,'gene')
                mut_annotations(i).gene2       = cdf.gene ;
            else
                mut_annotations(i).gene2 = '';
            end
            mut_annotations(i).locustag2 = cdf.locustag;
            mut_annotations(i).distance2 = cdf.loc1 - Positions(i,2);
            mut_annotations(i).protein2    = cdf.product ;
            if cdf.strand==0
                mut_annotations(i).distance2 = mut_annotations(i).distance2 * -1;
            end
        end
        mut_annotations(i).NonSyn = nan;
        mut_annotations(i).AA   = nan ;
    end
    
end

return
