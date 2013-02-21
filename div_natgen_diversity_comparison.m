function div_natgen_diversity_comparison


div_natgen_diversity_comparision

%Load
NGfixed=ismember(p,NGfixedp);
NGfixed(ancnt~=maNT(:,1))=0; %Removes positions that were called as fixed relative to isolates from Pt 1 but were really fixed in Pt 1
%NGfixed(ancnt~=maNT(:,1) & ancnt==maNT(:,comparisonpt) & maf(:,comparisonpt)>.95)=0; %Removes positions that were called as fixed relative to isolates from Pt 1 but were really fixed in Pt 1

Pfixed=fixedpos(:,comparisonpt)>0;


Ppolymorphic=goodpos(:,comparisonpt);
NGpolymorphic=ismember(p,NGpolymorphicp);
PpolymorphicT=goodpos(:,comparisonpt+1);%Later timepoint


%Fixed comparisons
disp([sum(NGfixed) sum(Pfixed)]) 
sum(NGfixed & Pfixed) %6, all that were fixed in NG!
sum(~NGfixed & Pfixed) %12, 12 new fixed mutations,
sum(NGpolymorphic & Pfixed) %... of which 4 were previously polymorphic and polymorphic only on account of the earliest isolate
%the remaining 8 are not shown as 'NGfixed' -- 7 are fixed in additional
%patients, 1 wasn't called in NG because there is no coverage in patient 1
%(deletion)
%so its a total of 7 specific mutations, which matches the supplemental
%figure


disp([sum(Ppolymorphic) sum(NGpolymorphic)]) 
sum(Ppolymorphic & NGpolymorphic) %=7, not many polymorphisms have hung around (out of 58 & 41)
sum(Ppolymorphic & PpolymorphicT & NGpolymorphic) %same as above-- no additional polymorphisms in later timepoint that were found earlier
sum(~Ppolymorphic & NGpolymorphic &~Pfixed) %24, Many polymorphic positions lost in population to time or detectability
sum(Ppolymorphic & ~NGpolymorphic) %=51, Many polymorphisms not previously found, at least one of which is not found due to thresholding

%%What are these 24 in NatGen but not in deep?
%polymorphicCallsNatGen(ismember(chrpos2index(polymorphicPositionsNatgen',ChrStarts),p(~Ppolymorphic & NGpolymorphic &~Pfixed)),:)
%%All were found in exactly one strain

%Not investigating the remaining 51 that are polymorphic in the new but not
%in the old


sum(NGfixed & Ppolymorphic) %=0, means nothing is currently polymorphic that was previously fixed, makes sense
sum(Pfixed & NGpolymorphic) %=4, means some mutations that were previously polymorphic have now fixed-- were all polymorphic due to one old isolate only



%NatGen figure suggests that there are more than 4 mutations seperating J-9-11 from others... are these not fixed also?
%polymorphicCallsNatGen shows 6 such lines (see below)
[~,NGpolyNTs]=ismember(polymorphicCallsNatGen(:,4:size(polymorphicCallsNatGen,2)-3),'ATCG');
justearly=zeros(size(NGpolymorphic));
for i=1:size(NGpolyNTs,1)
    line=NGpolyNTs(i,:);
    line(line==0)=mode(line);
    if numel(unique(line(1:end-1))==1) & unique(line(1:end-1))~=line(end)
        justearly(find(p==NGpolymorphicp(i),1))=1;
        disp(line)
    end
end
%what about the other two?
p(justearly  & NGpolymorphic & ~Pfixed)

%There are 4 mutations that fixed-
%They are 3 nonsynonymous mutations -- orbA(1606), gyrA (2180), and a beta
%lactamase (3694)
%There is a mutation in an intergenic region upstream of a
%methyltransferase 
