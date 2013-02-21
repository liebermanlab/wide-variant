function ntpairing = div_analyze_nearby_nts(p, goodpos, GenomeLength, SampleNames, SampleDirs, ChrStarts, combined_annotation_all, readlength)

NTs='ATCG';

ntpairing=struct;
k=1;
for i=2:size(goodpos,2)
    
    disp(SampleNames(i).Sample)
    goodp=p(goodpos(:,i)>0);
    nextp = [goodp(2:end); GenomeLength] - goodp ;
    candidatep=find(nextp<readlength);
    [~,lower]=ismember(goodp(candidatep),p);
    [~,upper]=ismember(goodp(candidatep+1),p);
    
    disp([lower upper])

    if numel(lower>0)
        disp(p(lower))
        disp(p(upper))
        c=div_pileup_region(p(lower), p(upper), SampleDirs(i), ChrStarts);
        
        for j=1:numel(lower)
            ntpairing(k).hits=squeeze(c(j,:,:));
            ntpairing(k).sample=SampleNames(i).Sample;
            ntpairing(k).pos1=p(lower(j));
            ntpairing(k).pos2=p(upper(j));
            ntpairing(k).index1=lower(j);
            ntpairing(k).index2=upper(j);
            ntpairing(k).annotation1=[combined_annotation_all(lower(j)).gene ' ' combined_annotation_all(lower(j)).annotation];
            ntpairing(k).annotation2=[combined_annotation_all(upper(j)).gene ' ' combined_annotation_all(upper(j)).annotation];
            ntpairing(k).ref1=find(NTs==combined_annotation_all(lower(j)).ref);
            ntpairing(k).ref2=find(NTs==combined_annotation_all(upper(j)).ref);
            k=k+1;
        end
    end
end

end