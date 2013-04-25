function div_cov_window(an, center, window, cwindows, names, samplen, showlegend)

global CONTROLSAMPLE


%initialize
strand=['g-','k-'];
scrsz = get(0,'ScreenSize');
figure(821);
set(821,'Position',[scrsz(3)/2 scrsz(4) scrsz(3)/2.5 scrsz(4)/4]);clf;hold on;
othersamples=1:numel(names); othersamples(CONTROLSAMPLE)=[];



%display name of genes and start/stop
if numel(an.locustag)>0
    protein='';
    for i=1:size(an.protein,1)
        protein=[protein ' ' an.protein(i,:)];
    end
    title([an.gene ':' an.locustag ' : ' protein '....' num2str(an.loc1) '-' num2str(an.loc2) ' Strand' num2str(an.strand+1)])
    plot([an.loc1,an.loc1],[0,1.10*max(cwindows(:))],strand(1+an.strand))
    plot([an.loc2,an.loc2],[0,1.10*max(cwindows(:))],strand(2-an.strand))
else
    title(['Intergenic::' num2str(an.distance1) '-' an.locustag1 ' , ' num2str(an.distance2) '-' an.locustag2])
    if abs(an.distance1) < window
        x1=center-abs(an.distance1);
        plot([x1,x1],[0,1.10*max(cwindows(:))],strand(1+(an.distance1>0)))
    end
    if abs(an.distance2) < window
        x1=center+abs(an.distance2);
        plot([x1,x1],[0,1.10*max(cwindows(:))],strand(1+(an.distance2>0)))
    end
end




%plot cov in window

if numel(names)>10
    set(gca, 'ColorOrder', colormap(gray));
end
l=(center-window):(center+window);
o = plot(l,cwindows(:,othersamples), '-', 'LineWidth', 1); %plot other samples
c = plot(l,cwindows(:,CONTROLSAMPLE),'k-', 'LineWidth',2); %plot control in black

legendlist=[c]; legendnames={names{CONTROLSAMPLE}};
if samplen > 0
    s = plot(l,cwindows(:,samplen),'r-', 'LineWidth',2); %plot clicked in red
    othersamples(othersamples==samplen)=[];
    legendlist(end+1)=s;
    legendnames{end+1}=names{samplen};
end

%Plot center
plot([center,center],[0,1.10*max(cwindows(:))],'k:')




%legend, etc
xlim([l(1) l(end)])
ylim([0 max(cwindows(:))])

ylabel('Coverage')


%legend
if showlegend >0
    for i=1:numel(othersamples)
        legendlist(end+1)=o(i);
        legendnames{end+1}=names{othersamples(i)};
    end
end
legend(legendlist, legendnames,'Location', 'BestOutside')


end

