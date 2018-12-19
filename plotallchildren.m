function plotallchildren(current,parentx,parenty1,parenty2)

global genotypes;
global parent;
global runningploty_byx;
global compartmentfreqs;
global sizescale;


directchildren=find(parent==current);

temp=compartmentfreqs; temp(isnan(temp))=0;
cfreqs=sum(temp,2);

alldescendents=find(sum(genotypes(:,find(genotypes(current,:))),2)==sum(genotypes(current,:)) & ~(1:numel(parent)==current)');
alldescendents_freq=sum(cfreqs(alldescendents)); 

if alldescendents_freq < (.5 * cfreqs(current));
    spacer = max(1,sizescale* (.5*(cfreqs(current)-alldescendents_freq))/numel(alldescendents));
else
    spacer=sizescale*.05;
end


for j=1:numel(directchildren) %[numel(directchildren):-1:1];
    c=directchildren(j);
    genotype=genotypes(c,:);
    runningploty_byx(sum(genotype))=runningploty_byx(sum(genotype))+spacer;
    [newchildy1, newchildy2] = makepatch(c,parentx+.5,parenty1,parenty2);
   % disp(newchildy1)
   % disp(newchildy2)
   % text(sum(genotype),runningploty,num2str(c)) %for troubleshooting
    if ~isempty(find(parent==c))
        plotallchildren(c,sum(genotype),newchildy1,newchildy2);
    end
end

if runningploty_byx(sum(genotypes(current,:))+1) < (parenty2+sizescale*.1);
    runningploty_byx(sum(genotypes(current,:))+1) = (parenty2+sizescale*.05);
end







