function div_maf_window(an, center, window, iwindow, dwindows, goodp, goodf, names, samplen, showlegend)



global CONTROLSAMPLE


%initialize
strand=['g-','k-'];
scrsz = get(0,'ScreenSize');
figure(820);
set(820,'Position',[scrsz(3)/20 scrsz(4) scrsz(3)/2.5 scrsz(4)/4]);clf;hold on;
legendlist=[];
legendnames={};
l=(center-window):(center+window);
othersamples=1:numel(names); othersamples(CONTROLSAMPLE)=[];

%display name of genes and start/stop
if numel(an.protein)>0
    if ~isfield(an,'locustag')
        an.locustag='';
    end
    protein='';
    for i=1:size(an.protein,1)
        protein=[protein ' ' an.protein(i,:)];
    end
    title([an.gene ':' an.locustag ' : ' protein '....' num2str(an.loc1) '-' num2str(an.loc2) ' Strand' num2str(an.strand+1)])
    plot([an.loc1,an.loc1],[0,1],strand(1+an.strand))
    plot([an.loc2,an.loc2],[0,1],strand(2-an.strand))
else
    title(['Intergenic::' num2str(an.distance1) '-' an.locustag1 ' , ' num2str(an.distance2) '-' an.locustag2])
    if abs(an.distance1) < window
        x1=center-abs(an.distance1);
        plot([x1,x1],[0,1],strand(1+(an.distance1>0)))
    end
    if abs(an.distance2) < window
        x1=center+abs(an.distance2);
        plot([x1,x1],[0,1],strand(1+(an.distance2>0)))
    end
end



%plot isolates, if isolate window
if (numel(iwindow)) >1
    
    %plot locations where major allele disagrees in different color --
    %indicated in iwindow by a negative allele frequency
    reverselocations=iwindow(2,:)<0;
    if sum(~reverselocations)>0
        iso=plot(iwindow(1,~reverselocations)-1,iwindow(2,~reverselocations),'vk','MarkerFaceColor', rgb('Orange'),'MarkerSize', 10);
        legendlist{end+1}={'isolates(S2)'};
    end
    if sum(reverselocations)>0
        iso=plot(iwindow(1,reverselocations)-1,-1*iwindow(2,reverselocations),'vk','MarkerFaceColor', rgb('LightGreen'),'MarkerSize', 10);
        legendlist{end+1}={'isolates(S2)- opposite major allele'};
    end
    legendlist(end+1)=iso;

end


%plot major allele frequency in window
if numel(names)>10
    set(gca, 'ColorOrder', colormap(gray));
end

o = plot(l,dwindows(:,othersamples), '.', 'MarkerSize', 15); %plot other samples
c = plot(l,dwindows(:,CONTROLSAMPLE),'k.', 'MarkerSize', 15); %plot control in black

legendlist(end+1)=c; legendnames{end+1}=names{CONTROLSAMPLE};
if samplen > 0
    s = plot(l,dwindows(:,samplen),'r.', 'MarkerSize', 15); %plot sample in red
    legendlist(end+1)=s;
    legendnames{end+1}=names{samplen};
    othersamples(othersamples==samplen)=[];
end


%display diamonds at good positions -give solid backdrop and replot
%original dot
plot(goodp, goodf,'kd', 'MarkerSize', 10, 'MarkerFaceColor', rgb('Silver'));
plot(goodp, goodf,'.', 'MarkerSize', 15)
if samplen~=CONTROLSAMPLE
    plot(goodp, goodf(:,samplen),'kd', 'MarkerSize', 10,'MarkerFaceColor', rgb('Pink'));
    plot(goodp, goodf(:,samplen),'r.', 'MarkerSize', 15)
end
    

%Plot center
plot([center,center],[0,1.1],'k:')





%legend
if showlegend >0
    for i=1:numel(othersamples)
        legendlist(end+1)=o(i);
        legendnames{end+1}=names{othersamples(i)};
    end
end
legend(legendlist, legendnames,'Location', 'BestOutside')



%label
axis([l(1) l(end) .5 1.01])
ylabel('Major allele frequency')

end

