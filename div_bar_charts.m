function div_bar_charts(c, samplen, names, params)



%params is an optional argument to display current threshold settings

name=names{samplen};

scrsz = get(0,'ScreenSize');

figure(660);clf;
%set(660,'Position',[scrsz(3)*2/3 scrsz(4)/20 scrsz(3)/3 scrsz(4)/2]);clf;hold on;


%Plot number of calls
subplot(4,1,1); hold on;
title(['Number of calls ' ])%.... ' name ' p(no strand bias) = ' num2str(c(end-6,samplen))]);
a=squeeze(c(1:8,:));
bar(reshape([a; nan(4,size(a,2))],4,[])','stacked')
legend('A','T','C','G', 'Location', 'BestOutside')
ylabel('Number of reads')
set(gca,'Xtick',2:3:(3*numel(names)-1))
set(gca,'Xticklabel',names)
xlim([0 (3*numel(names)+3)])
xticklabel_rotate;

%Plot call quality
subplot(4,1,2); hold on;
title(['Average call quality ' ])%.... ' name ' p = ' num2str(c(end-5,samplen))]);
a=squeeze(c(9:16,:));
if nargin>3
    plot([0 3*size(a,2)], [params.min_bq params.min_bq], 'k:')
end
bar(reshape([a; nan(4,size(a,2))],4,[])','grouped', 'LineStyle', 'none')
legend('Aq','Tq','Cq','Gq', 'Location', 'BestOutside')
ylabel('Average Base Quality')
set(gca,'Xtick',2:3:(3*numel(names)-1))
set(gca,'Xticklabel',names)
xlim([0 (3*numel(names)+3)])
xticklabel_rotate;


%Plot mapping quality
subplot(4,1,3); hold on;
title(['Average mapping quality ' ])% .... ' name ' p = ' num2str(c(end-4,samplen))]);
a=squeeze(c(17:24,:));
if nargin>3
    plot([0 3*size(a,2)], [params.min_mq params.min_mq], 'k:')
end
bar(reshape([a; nan(4,size(a,2))],4,[])','grouped','LineStyle', 'none')
legend('Am','Tm','Cm','Gm', 'Location', 'BestOutside')
ylabel('Average Mapping Quality')
set(gca,'Xtick',2:3:(3*numel(names)-1))
set(gca,'Xticklabel',names)
xlim([0 (3*numel(names)+3)])
xticklabel_rotate;

%Plot tail distance f
subplot(4,1,4); hold on;
title(['Average tail distance  ' ]) % (each strand).... ' name ' fp = ' num2str(c(end-3,samplen))  ' rp = ' num2str(c(end-2,samplen))]);
a=squeeze(c(25:32,:));
if nargin>3
    plot([0 3*size(a,2)], [params.min_td params.min_td], 'k:')
    plot([0 3*size(a,2)], [params.max_td params.max_td], 'k:')
end
bar(reshape([a; nan(4,size(a,2))],4,[])','grouped', 'LineStyle', 'none')
legend('Atd','Ttd','Ctd','Gtd', 'Location', 'BestOutside')
ylabel('Average Tail Distance')
set(gca,'Xtick',2:3:(3*numel(names)-1))
set(gca,'Xticklabel',names)
xlim([0 (3*numel(names)+3)])
xticklabel_rotate;


end

