function div_maf_window_metagenomics(marker, pos, windows, goodp, goodf, names)


%initialize
strand=['g-','k-'];
scrsz = get(0,'ScreenSize');
figure(820);

set(820,'Position',[1 scrsz(4) scrsz(3) scrsz(4)/3]);clf;hold on;



%plot maf in window
plot(windows, '.', 'MarkerSize', 15) %plot all samples

%display diamond at good positions
plot(goodp, goodf,'kd', 'MarkerSize', 10);

%indicate center of window
plot([pos,pos],[0,1],'r')

title(marker)

%legend, etc
axis([l(1) l(end) .5 1])
ylabel('MAF')
legend(names);


end

