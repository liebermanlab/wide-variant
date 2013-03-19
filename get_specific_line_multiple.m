function get_specific_line_multiple(filename,SampleName, ScafNames,mut_list)


%Modification to allow for search for all lines in strain.vcf that match
%that position. Only record info for best position.

Tab=9 ;


fid = fopen(filename,'r') ;
lins = {} ;
indx = [] ;
pos_list = mut_list(:,1)*1e8 + mut_list(:,2) ;
[pos_list,ind] = sort(pos_list) ;

i=0 ;
ipos=1 ;
topscore=0;
prevpos=0;
y=0 ; %position in linds, lindx
tic
while ~feof(fid) && ~isempty(ipos)
    i=i+1 ;
    if ~mod(i,50000), fprintf(1,'.'); end
    s = fgetl(fid) ;
    f = find(s==Tab,6) ;
    if length(f)>=2
        raw_pos = str2num(s(f(1)+1:f(2)-1)) ;
        if isempty(raw_pos), raw_pos=0 ; end
        scf = s(1:f(1)-1) ;
        iscf = find(strcmp(scf,ScafNames)) ;
        pos = iscf*1e8 + raw_pos ;
        if pos>=pos_list(ipos)
            if any(pos==pos_list)
                score=str2num(s(f(5)+1:f(6)-1));  %right now this uses the score column, might be best to use fq column... though this is used in downstream processing
                if prevpos~=pos
                    y=y+1;
                    topscore=0;
                    ipos = find(pos_list>prevpos,1) ;
                    prevpos=pos ;
                end
                if score >= topscore &  ~(s(f(4):f(4)+1)=='.' &  strcmp(s(f(6)+3:f(6)+7),'INDEL')) %%write only if new or better than previous guess.  some cases where samtools shows a weakly supported indel as a variant without putting it into the column
                    lins{y} = s ;
                    indx(y) = find(pos_list==pos,1) ;
                    topscore=score;
                end
            end
        end
    end
end
fprintf(1,' %6.0f \n',toc) ;
indx=ind(indx) ;

save(['vcfinfo_' char(SampleName)], 'lins', 'indx')

end

