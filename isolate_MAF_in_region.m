function [imaf, iminor, imajorNT, imaf_window]=isolate_MAF_in_region(dpos,ipos,maf, minor, majorNT, cstarts, wl)

%for both spos and ipos, each column is a position
%first row is column, second row is position

Npos=numel(dpos);

imaf=ones(Npos,1); %maf of isolate at position, 1 otherwise
iminor=zeros(Npos,1); %minor of isolate at position, 1 otherwise
imajorNT=zeros(Npos,1); %major NT of isolates at position 1, 0 otherwise

[inter,loc]=ismember(dpos,ipos);
imaf(inter)=maf(loc(loc>0));
iminor(inter)=minor(loc(loc>0));
imajorNT(inter)=majorNT(loc(loc>0));

imaf_window=cell(Npos,1); %in nearby region

for p=1:Npos
    region=dpos(p)-wl:dpos(p)+wl;
    [~,loc]=ismember(region,ipos);
    if ~isempty(loc(loc>0))
        temp=zeros(2,numel(loc(loc>0)));
        positions=p2chrpos(ipos(loc(loc>0)),cstarts);
        temp(1,:)= positions(:,2);
        temp(2,:)=maf(loc(loc>0)); %maf
        imaf_window(p)={temp};
    end
end

