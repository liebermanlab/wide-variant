function div_data_hists(d, nontest_reqs, p, savename)

scrsz = get(0,'ScreenSize');
num_nontest_reqs=max(nontest_reqs);

%unpack parameters
minorfreqthreshold=p.minorfreqthreshold;
maxreads_perstrand=p.maxreads_perstrand;
minreads_perstrand=p.minreads_perstrand;
min_bq=p.min_bq;
min_mq=p.min_mq;
max_td=p.max_td;
min_td=p.min_td;
max_sbp=p.max_sbp;
max_bqp=p.max_bqp;
%max_mqp=p.max_mqp;
max_tdp=p.max_tdp;
max_id=p.max_percent_indels;

%get p values
SB=d(end-6,:);
BQ=d(end-5,:);
MQ=d(end-4,:);
TDF=d(end-3,:);
TDR=d(end-2,:);


%Before plotting, replace extremely high values of p values
SB(SB>40)=40;
BQ(BQ>40)=40;
MQ(MQ>40)=40;
TDF(TDF>40)=40;
TDR(TDR>40)=40;


[maf, majorNT, minorNT] = div_major_allele_freq(d);
Npos=size(d,2);



%Create histograms of d
h1=figure('Position',[1 1 scrsz(3)/1.5 scrsz(4)]);clf;hold on; clf;

subplot(4,4,1); hold on;
hist(sum(d(1:8,:).*d(9:16,:))./sum(d(1:8,:)),50);
n=hist(sum(d(1:8,:).*d(9:16,:))./sum(d(1:8,:)),50);
xlabel('Mean base quality')
ylabel('Genomic positions')
plot([min_bq min_bq], [0 max(n)], 'g-','LineWidth',2)
axis([0 50 0 max(n)])


subplot(4,4,2); hold on;
hist(sum(d(1:8,:).*d(17:24,:))./sum(d(1:8,:)),50);
n=hist(sum(d(1:8,:).*d(17:24,:))./sum(d(1:8,:)),50);
xlabel('Mean mapping quality')
ylabel('Genomic positions')
plot([min_mq min_mq], [0 max(n)], 'g-','LineWidth',2)
axis([0 50 0 max(n)])


subplot(4,4,3); hold on;
hist(sum(d(1:4,:).*d(25:28,:))./sum(d(1:4,:)),50);
n=hist(sum(d(1:4,:).*d(25:28,:))./sum(d(1:4,:)),50);
xlabel('Mean tail distance -forward strand')
ylabel('Genomic positions')
plot([min_td min_td], [0 max(n)], 'g-','LineWidth',2)
plot([max_td max_td], [0 max(n)], 'g-','LineWidth',2)
axis([0 50 0 max(n)])

subplot(4,4,4); hold on;
hist(sum(d(5:8,:).*d(29:32,:))./sum(d(5:8,:)),50);
n=hist(sum(d(5:8,:).*d(29:32,:))./sum(d(5:8,:)),50);
xlabel('Mean tail distance -reverse strand')
ylabel('Genomic positions')
plot([min_td min_td], [0 max(n)], 'g-','LineWidth',2)
plot([max_td max_td], [0 max(n)], 'g-','LineWidth',2)
axis([0 50 0 max(n)])


subplot(4,4,5);  hold on;
hist(double(SB(nontest_reqs==num_nontest_reqs)),50);
n=hist(double(SB(nontest_reqs==num_nontest_reqs)),50);
xlabel('-log p value, strand bias')
ylabel('Positions meeting nontest thresholds')
plot([max_sbp max_sbp], [0 max(n)], 'g-','LineWidth',2)
axis([0 40 0 max(n)])

subplot(4,4,6); hold on;
hist(double(BQ(nontest_reqs==num_nontest_reqs)),40);
n=hist(double(BQ(nontest_reqs==num_nontest_reqs)),40);
xlabel('-log p value, base quality')
ylabel('Positions meeting nontest thresholds')
plot([max_bqp max_bqp], [0 max(n)], 'g-','LineWidth',2)
axis([0 40 0 max(n)])

subplot(4,4,7); hold on;
hist(double(MQ(nontest_reqs==num_nontest_reqs)),40);
n=hist(double(MQ(nontest_reqs==num_nontest_reqs)),40);
xlabel('-log p value, mapping quality')
ylabel('Positions meeting nontest thresholds')
%plot([max_mqp max_mqp], [0 max(n)], 'g-','LineWidth',2)
axis([0 40 0 max(n)])

subplot(4,4,8); hold on;
hist(double(TDF(nontest_reqs==num_nontest_reqs)),40);
n=hist(double(TDF(nontest_reqs==num_nontest_reqs)),40);
xlabel('-log p value, tail distance F')
ylabel('Positions meeting nontest thresholds')
plot([max_tdp max_tdp], [0 max(n)], 'g-','LineWidth',2)
axis([0 40 0 max(n)])


subplot(4,4,9); hold on;
hist(double(TDR(nontest_reqs==num_nontest_reqs)),40);
n=hist(double(TDR(nontest_reqs==num_nontest_reqs)),40);
xlabel('-log p value, tail distance R')
ylabel('Positions meeting nontest thresholds')
plot([max_tdp max_tdp], [0 max(n)], 'g-','LineWidth',2)
axis([0 40 0 max(n)])


subplot(4,4,10); hold on;
hist(sum(d(1:4,:)),240);
n=hist(sum(d(1:4,:)),240);
xlabel('Forward strand coverage')
ylabel('Positions')
plot([maxreads_perstrand maxreads_perstrand], [0 max(n)], 'g-','LineWidth',2)
plot([minreads_perstrand minreads_perstrand], [0 max(n)], 'g-','LineWidth',2)
axis([0 1200 0 5*10^5])


subplot(4,4,11);  hold on;
hist(sum(d(5:8,:)),240);
n=hist(sum(d(5:8,:)),240);
xlabel('Reverse strand coverage')
ylabel('Positions')
plot([maxreads_perstrand maxreads_perstrand], [0 max(n)], 'g-','LineWidth',2)
plot([minreads_perstrand minreads_perstrand], [0 max(n)], 'g-','LineWidth',2)
axis([0 1200 0 5*10^5])


subplot(4,4,12);  hold on;
percent_indels=double(d(end,:))./(sum(d(1:8,:))+double(d(end,:)));
[x,n]=hist(percent_indels,50);
c=log10(x); bar(n,c);
xlabel('Percent reads supporting indels')
ylabel('Log10 Positions')
plot([max_id max_id], [0 max(c)], 'g-','LineWidth',2)
axis([0 1 0 max(c)])

subplot(4,4,13);  hold on;
minorfreq=double(d(sub2ind(size(d),minorNT,1:Npos))+d(sub2ind(size(d),minorNT+4,1:Npos)))./sum(d(1:8,:));
[x,n]=hist(minorfreq,50);
c=log10(x); bar(n,c);
xlabel('Minor allele frequency')
ylabel('Log10 Positions')
plot([minorfreqthreshold minorfreqthreshold], [0 max(c)], 'g-','LineWidth',2)
axis([0 .5 0 max(c)])

print(h1, savename)

end