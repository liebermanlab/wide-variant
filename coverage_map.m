function coverage_map(dir)
genome_length = 2821361;
raw = zeros(genome_length,1);

window = 500;

fid = fopen([dir '/strain.pileup'],'r');

data = textscan(fid,'%s %d %c %d %*[^\n]','delimiter','\t');

pos = data{2};
for i = 1:length(pos)
    raw(pos(i)) = data{4}(i);
end

smoothed = smooth(raw,window);

save([dir '/coverage.mat'],'raw','smoothed');

f=figure;
plot(raw,'b');
xlabel('Position');
ylabel('Read Depth');

print(f,[dir '/coverage_map'],'-dtiff');
close(f);

h=figure;
hist(raw,0:25:500);
xlabel('Read Depth');
ylabel('Base Pairs');

print(h,[dir '/coverage_hist'],'-dtiff');
close(h);

end