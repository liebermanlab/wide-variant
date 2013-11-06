function lung_get_all_clonality(mut_freq, SampleNames, isolates_per_site)
    
    if nargin < 3
        isolates_per_site = 24; 
    end

    % divide isolate numbers per site 
    strbit = '-'; 
    SiteNames = get_unique_site_names(SampleNames); 
    Sites = cell(length(SiteNames),1); 
    for n = 1:length(SiteNames)
        Sites{n} = (isolates_per_site*(n-1)+1):isolates_per_site*n;
    end
    
    
    % INITIALIZE DATA STRUCTURE 
    
    % number of unique clones wt/bt sites
    clonality_matrix = zeros(length(Sites)); 
    % number of SHARED clones bt sites
    shared_clones_matrix = zeros(length(Sites)); 
    % per isolate, what "type" is it within site
    all_isolate_types = cell(length(Sites),2); 
    % number of isolates per clonal type within site
    unique_isolates = cell(length(Sites),1);
    
    % WITHIN SITE 
    for i = 1:length(Sites)
        site_isolates = mut_freq(:,Sites{i}); 
        [isolate_type, cnts_per_type, site_unique_genotypes] = lung_calc_clonal_isolates(site_isolates); 
        % stats on number of unique types
        all_isolate_types{i,1} = isolate_type;
        all_isolate_types{i,2} = cnts_per_type;
        % # of unique clones
        clonality_matrix(i,i) = length(cnts_per_type); 
        % store unique genotypes per site
        unique_isolates{i} = site_unique_genotypes;
    end

    % BETWEEN SITES
    for si = 1:length(Sites)
        for sj = 1:length(Sites)
            genotype_si = unique_isolates{si}; 
            genotype_sj = unique_isolates{sj}; 
            % get # of unique clones 
            [i_type, ~, ~] = lung_calc_clonal_isolates([genotype_si, genotype_sj]); 
            num_unique_clones = length(unique(i_type));
            
            shared_clone_cnt = size(genotype_si,2) + size(genotype_sj,2) - num_unique_clones;
            
            if si ~= sj
                clonality_matrix(si,sj) = num_unique_clones; 
                clonality_matrix(sj,si) = num_unique_clones; 
                shared_clones_matrix(si,sj) = shared_clone_cnt; 
                shared_clones_matrix(sj,si) = shared_clone_cnt; 
            end
        end
    end
    
    % PLOT # OF CLONAL ISOLATES PER SITE
    plot_clonality_matrix(tril(clonality_matrix), 'Number of Clonal Groups'); 
    
    % KMEANS CLUSTERING
    cluster_size = 3; 
    % clonality_matrix is not the best to cluster on 
    % interested in CO-OCCURENCE of sites with each other 
    clusterindex = kmeans(clonality_matrix, cluster_size); 
    for i = 1:cluster_size
        fprintf('\nCluster group %i', i); 
        SiteNames(clusterindex==i)
    end
    
    % PLOT # OF SHARED CLONES BETWEEN SITES 
    plot_clonality_matrix(tril(shared_clones_matrix), 'Number of Shared Isolates'); 
    

    % ___ Subfunctions ___ %
    function uniquekeys = get_unique_site_names(snames)
        allkeys = {}; 
        for u = 1:length(snames)
            sn = snames{u}; 
            keyloc = strfind(sn,strbit);
            skey = sn(1:keyloc-1); 
            allkeys{end+1} = skey;
        end
        uniquekeys = unique(allkeys, 'stable'); % prevents sorting 
    end

    function plot_clonality_matrix(clonal_counts, titlename)
        figure, imagesc(clonal_counts); 
        title(titlename, 'FontSize', 24, 'FontWeight', 'bold'); 
        colormap(flipud(gray)); 

        textStrings = num2str(clonal_counts(:), '%i'); 
        textStrings = strtrim(cellstr(textStrings)); 
        [x,y] = meshgrid(1:length(SiteNames)); 
        hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
                        'HorizontalAlignment','center', ...
                        'FontSize', 16);
        midValue = mean(get(gca,'CLim'));  %# Get the middle value of the color range
        textColors = repmat(clonal_counts(:) > midValue,1,3); 
        set(hStrings,{'Color'},num2cell(textColors,2));
        set(gca, 'XTick', 1:length(SiteNames), ...
                'YTick', 1:length(SiteNames), ...
                'XTickLabel', SiteNames, ...
                'YTickLabel', SiteNames, ...
                'FontSize', 20)
    end

end