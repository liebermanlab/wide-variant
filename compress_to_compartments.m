function compartmentm = compress_to_compartments(matrix,compartmentlist)

global compartments;

compartmentm=zeros(size(matrix,1),numel(compartments));

%disp(size(matrix))
%disp(size(compartmentm))
for c=1:numel(compartments)
    if ismember(c,compartmentlist);
        numsamplesc=sum(compartmentlist==c);
     %   disp(size(sum(matrix(:,compartmentlist==c),2)/numsamplesc))
      %  disp(size(compartmentm(c,:)))
        compartmentm(:,c)=sum(matrix(:,compartmentlist==c),2)/numsamplesc;
    end
end
        
