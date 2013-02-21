function normalized_specificity = div_analyze_specificity_using_isolates(iAF,pAF, cnts, params, covthresholds, sample, figuren)

%December 2012 TDL
%Specificity returned is using specified parameters. Graphs explore allele
%frequency cutoff



%get data passing all other thresholds without respect to allele frequency
no_freq_params=params;
no_freq_params.minorfreqthreshold=.00;
testpos=div_test_thresholds(cnts, no_freq_params, covthresholds);

%observed specificity by rank
pooledf=pAF(testpos(:,sample)>0,sample);
isolatesf=iAF(testpos(:,sample)>0);
[pooledf, sortedpositions]=sort(pooledf,1,'descend');
isolatesf=isolatesf(sortedpositions);
ospecificity=cumsum(isolatesf>0)./cumsum(ones(size(isolatesf)));

%expected specificity by rank
prob_in_isolates=1.-((1-pooledf).^29);
especificity=cumsum(prob_in_isolates)./cumsum(ones(size(isolatesf)));


freqcutoff=find(pooledf <= params.minorfreqthreshold, 1)-1;
% disp(pooledf)
% %disp(params.minorfreqthreshold)
if isempty(freqcutoff)
    normalized_specificity=0;
else
   % fprintf(num2str(freqcutoff+1))
    %fprintf('blahadfa')
    if figuren > 0
        %plot in terms of absolute
        figure(figuren); clf; hold on;
        %background
        freqcutoffs=[.15 .08 .05 .03 .02 .015]; lowerx=0;
        for i=1:numel(freqcutoffs)
            upperx=find(pooledf <= freqcutoffs(i), 1)-1;
            p1=patch([lowerx; upperx; upperx; lowerx], [0.01; 0.01; numel(prob_in_isolates); numel(prob_in_isolates)],  repmat(1-5*freqcutoffs(i),1,3));
            lowerx=upperx; set(p1,'EdgeColor','None','LineWidth',.01)
        end
        plot(1:numel(prob_in_isolates), 1:numel(prob_in_isolates), 'k:','LineWidth',1.1)
        plot(cumsum(prob_in_isolates), 'k', 'LineWidth',4)
        plot(cumsum(isolatesf>0),'m', 'LineWidth',4)
        xlabel('Rank polymorphic position in pooled P1')
        ylabel('Number of mutations also found in 29 isolates from P1')
        plot([freqcutoff freqcutoff], [0 max(cumsum(prob_in_isolates))], 'r--','LineWidth',3) %bars to indicate cutoff
        axis([0 upperx 0 find(pooledf <= params.minorfreqthreshold -.015, 1)]);
       % set(gca,'XScale','log')
        
        
        %plot normalized and smoothed with a window of 5
        figure(100+ figuren); clf; hold on;
        freqcutoffs=[.15 .08 .05 .03 .02 .015]; lowerx=0;
        for i=1:numel(freqcutoffs)
            upperx=find(pooledf <= freqcutoffs(i), 1)-1;
            p1=patch([lowerx; upperx; upperx; lowerx], [0.001; 0.001; 2; 2],  repmat(1-5*freqcutoffs(i),1,3));
            lowerx=upperx; set(p1,'EdgeColor','None','LineWidth',.01)
        end
        plot(smooth(ospecificity./especificity,5),'m', 'LineWidth',4)
        xlabel('Rank polymoprhic position found in pooled P1')
        ylabel('Specificity observed / Specificiy expected')
        plot([freqcutoff freqcutoff], [0 1.5],'r--','LineWidth',3) %bars to indicate cutoff
        axis([0 upperx 0 1.02]);
        
    end
    
    x=ospecificity./especificity;
    % disp(ospecificity)
    % disp('blah')
    % disp(['O/E, O, E'])
    % disp([x(freqcutoff) ospecificity(freqcutoff) especificity(freqcutoff)] )
    %disp(x)
    %disp(freqcutoff)
    normalized_specificity=x(freqcutoff);
    
end
