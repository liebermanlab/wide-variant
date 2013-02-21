% coverage_analysis.m
%
% Seungsoo Kim
% January 3, 2013
%
% analysis of coverage - determines breakpoints of amplifications/deletions

isolates=read_isolate_table;

if ~exist('coverage_processed.mat')
    load('./WT/trim_end25_20/bowtie2/coverage.mat');
    wt = raw/mean(raw);
    samples = zeros(2821361,91);
    deletions = zeros(2821361,91);

    for i = 1:length(isolates)
        s=isolates(i).Sample;
        load(['./' s '/trim_end25_20/bowtie2/coverage.mat']);
        samples(:,i) = raw/mean(raw)-wt;
        deletions(:,i) = raw == 0;
    end

    save('coverage_processed.mat','samples','wt','deletions');
else
    figure;hold on;
    for i = 1:91
        plot(find(samples(:,i)>1.5),(92-i)*ones(length(find(samples(:,i)>1.5)),1),'b.');
        plot(find(samples(:,i)<-1.5),(92-i)*ones(length(find(samples(:,i)<-1.5)),1),'r.');
        xlim([1 2821361]);
    end
end