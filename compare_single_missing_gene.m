function pairwise_std = compare_single_missing_gene(genenum, gene_coverage, gene_error, genes_nr, ...
                                                        allcoverage, SampleNames, allmodes)

    cov = gene_coverage(genenum,:); 
    err = gene_error(genenum,:); 
    pairwise_std = zeros(length(cov)); 
    
    for i = 1:length(cov)
        for j = 1:length(cov) 
            pairwise_std(i,j) = (cov(i)-cov(j))./sqrt(err(i).^2+err(j).^2); 
        end
    end
    
    [maxval,ind] = max(pairwise_std(:)); 
    [colony_i, colony_j] = ind2sub(length(SampleNames),ind);

    
%     if maxval > 0 
%         
%         genepos = genes_nr(genenum).indices;
%         if genepos(1)>genepos(2)
%             genepos = flipdim(genepos,2);
%         end
%         
%         % plot on fixed figure handle
%         scrsz = get(0,'ScreenSize'); 
%         figure(101); clf; hold on;
%         set(101,'Position',[scrsz(3)*1/3 scrsz(4) scrsz(3)/3 scrsz(4)/4]);
%         
%         if genepos(1) <= 100
%             xpos = (genepos(1)):genepos(2)+100; 
%             front_offset = 0;
%             end_offset = 100; 
%         else
%             xpos = (genepos(1)-100):(genepos(2)+100); 
%             front_offset = 100; 
%             end_offset = 100; 
%         end
%         plot(allcoverage(colony_i, xpos)'./allmodes(colony_i), 'b', 'LineWidth', 1.5); 
%         plot(allcoverage(colony_j, xpos)'./allmodes(colony_j), 'k', 'LineWidth', 1.5); 
%         legend(SampleNames{colony_i}, SampleNames{colony_j}, ...
%                 'FontSize', 16, 'FontWeight', 'bold', 'Location', 'NorthWest'); 
%         
%         % plot params
%         genelen = genepos(2)-genepos(1); 
%         
%         % plot gene limits 
%         plot([front_offset front_offset], [0 1000], 'r'); 
%         plot([genelen+front_offset genelen+front_offset], [0 1000], 'g'); 
%         set(gca, 'XTick', 1:500:genelen+front_offset+end_offset, ...
%                 'XTickLabel', genepos(1)-front_offset:500:genepos(2)+end_offset); 
%         xlim([0 genelen+front_offset+end_offset]); 
%         ylim([0 5]);
%         title(sprintf('Gene %i, Std %0.2f, %s', genenum, maxval, genes_nr(genenum).product), ...
%             'FontWeight', 'bold', 'FontSize', 16); 
%         xlabel('Genome Position', 'FontSize', 16, 'FontWeight', 'bold'); 
%         ylabel('# Reads/Mode', 'FontSize', 16, 'FontWeight', 'bold'); 
%         
%     end
    
    
end