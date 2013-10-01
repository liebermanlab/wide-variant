%% manual look for P09 trial 
P09_SP1 = ~cellfun(@isempty, strfind(samplenames_threshold, '09SP1-')); 
P09_SP4 = ~cellfun(@isempty, strfind(samplenames_threshold, '09SP4-')); 
P09_ST1 = ~cellfun(@isempty, strfind(samplenames_threshold, '09ST1-')); 

% just the purity
figure, imagesc(missing_isolates(:,1:3)); 
colormap('bone'); 
set(gca, 'XTick', 1:3, ...
        'XTickLabel', {'SP1-Pur', 'SP4-Pur1', 'SP4-Pur2'}, ...
        'YTick', [], ...
        'TickLength', [0 0], ... 
        'FontSize', 12); 
title('Patient 09 Purity', 'FontSize', 24, 'FontWeight', 'bold'); 
    

figure, imagesc(missing_isolates(:,P09_SP1)); colormap('bone'); 
title('P09 SPUTUM 1', 'FontSize', 24, 'FontWeight', 'bold'); 
set(gca, 'YTick', [], 'XTick', []); 
xlabel('Isolates', 'FontSize', 20, 'FontWeight', 'bold'); 

figure, imagesc(missing_isolates(:,P09_SP4)); colormap('bone'); 
title('P09 SPUTUM 4', 'FontSize', 24, 'FontWeight', 'bold'); 
set(gca, 'YTick', [], 'XTick', []); 
xlabel('Isolates', 'FontSize', 20, 'FontWeight', 'bold'); 

figure, imagesc(missing_isolates(:,P09_ST1)); colormap('bone'); 
title('P09 STOOL 1', 'FontSize', 24, 'FontWeight', 'bold'); 
set(gca, 'YTick', [], 'XTick', []); 
xlabel('Isolates', 'FontSize', 20, 'FontWeight', 'bold'); 

%% clustergram
cg1 = clustergram(missing_isolates(:,P09_SP1), 'Colormap', redbluecmap); 
t1 = addTitle(cg1, 'P09 SPUTUM 1'); 
set(t1, 'FontSize', 24', 'FontWeight', 'bold'); 

cg2 = clustergram(missing_isolates(:,P09_SP4), 'Colormap', redbluecmap); 
t2 = addTitle(cg2, 'P09 SPUTUM 4'); 
set(t2, 'FontSize', 24', 'FontWeight', 'bold'); 

cg3 = clustergram(missing_isolates(:,P09_ST1), 'Colormap', redbluecmap); 
t3 = addTitle(cg3, 'P09 STOOL 1'); 
set(t3, 'FontSize', 24', 'FontWeight', 'bold'); 