function div_plot_dnds(fignum, bins, cutoff, N, S, samplename, expected, upperexpected, lowerexpected)


figure(fignum);  clf;
scrsz = get(0,'ScreenSize');

set(fignum,'Position',[scrsz(3)/8 scrsz(4) scrsz(3)/4 scrsz(4)/2]);clf;hold on;


[~, alldNdS, ~] = div_calc_dnds_ci(N,S, expected);


nummuts=N+S;

i=find(nummuts>1,1,'first');
nummuts=nummuts(i:end);
N=N(i:end);
S=S(i:end);
bins=bins(i:end);
alldNdS=alldNdS(i:end);

subplot(2,1,1);  hold on;
plot([.5 1], [1 1], 'k-') %expected
plot([bins; bins],[lowerexpected(nummuts)./expected; upperexpected(nummuts)./expected], 'k-', 'LineWidth', 3)
plot(bins,alldNdS, 'rd', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
set(gca,'YScale','log')
set(gca,'YTick', [0.5 1 2 3 4])
axis([-Inf Inf 0.5 Inf])
xlabel('Major allele frequency cutoff (positions with MAF <x considered)')
ylabel('dN/dS (normalized to expectation)')





i=find(bins>=cutoff,1);
dNdS_diverse = div_calc_dnds(N(i),S(i), expected);
dNdS_less_diverse = div_calc_dnds(N(end)-N(i),S(end)-S(i), expected);

subplot(2,1,2);  hold on;
plot(1,dNdS_diverse, 'rd', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
plot(2,dNdS_less_diverse, 'bd', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
plot([1 1],[upperexpected(N(i)+S(i))./expected lowerexpected(N(i)+S(i))./expected]', 'k-', 'LineWidth', 3)
plot([2 2],[upperexpected(N(end)+S(end)-N(i)-S(i))./expected lowerexpected(N(end)+S(end)-N(i)-S(i))./expected]', 'k-', 'LineWidth', 3)
plot([0 3], [1 1], 'k-') %expected
axis([-Inf Inf 0.5 3])
set(gca,'YScale','log')
set(gca,'YTick', [0.5 1 2 3 4])
legend({['diverse: MAF <' num2str(cutoff)], ['less diverse: MAF >' num2str(cutoff)], 'expected 95% CI',}, 'Location', 'SouthWest')
ylabel('dN/dS (normalized to expectation)')
title(samplename)

end