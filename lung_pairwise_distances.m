function lung_pairwise_distances(mut_freq, SampleNames, isolates_per_site)
    % USES PDIST 

    if nargin < 3
        isolates_per_site = 24; 
    end
    
    % plot param
    max_dist = 60; 

    % divide isolate numbers per site 
    strbit = '-'; 
    SiteNames = get_unique_site_names(SampleNames); 
    Sites = cell(length(SiteNames),1); 
    for n = 1:length(SiteNames)
        Sites{n} = (isolates_per_site*(n-1)+1):(isolates_per_site*n);
    end

    % plot pairwise distance heatmap within each site (lower triangle!) 
    D = {}; 
    figure; 
    for i=1:length(SiteNames)
        dist = tril(squareform(pdist(mut_freq(:,Sites{i})', 'cityblock')));
        D{end+1} = dist;
        ax(i) = subplot(3,5,i); 
        imagesc(dist, [0 max_dist]); 
        title(SiteNames{i}, 'FontSize', 20, 'FontWeight', 'bold'); 
        set(gca, 'XTick', [], 'YTick', []); 
    end
   
    h=colorbar;
    set(h, 'Position', [.95 .11 .02 .7], 'Clim', [0 max_dist], 'YTick', [0 20 40 60], ...
            'FontSize', 20)
    for s=1:length(Sites)
        pos=get(ax(s), 'Position');
        set(ax(s), 'Position', [pos(1) pos(2) 0.9*pos(3) pos(4)]);
    end

    function uniquekeys = get_unique_site_names(snames)
        allkeys = {}; 
        for u = 1:length(snames)
            sn = snames{u}; 
            keyloc = strfind(sn,strbit);
            skey = sn(1:keyloc-1);
            allkeys{end+1} = skey;
        end
        uniquekeys = unique(allkeys, 'stable');
    end

end

