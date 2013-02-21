function sensitivity = div_analyze_sensitivity_using_isolates(iAF,pAF, cnts, params, covthresholds, sample, figuren)

%December 2012 TDL 
%Sensitivity returned is using specified parameters. Graphs explore allele
%isolate frequency cutoff



%get data passing all other thresholds without respect to allele frequency
no_freq_params=params;
no_freq_params.minorfreqthreshold=.00;
testpos=div_test_thresholds(cnts, no_freq_params, covthresholds);
pAF(testpos==0)=0; pooledf=pAF(:,sample);



%observed specificity by rank
isolatesf=iAF(iAF>0 & iAF <1);
pooledf=pooledf(iAF>0 & iAF <1);
[isolatesf, sortedpositions]=sort(isolatesf,1,'descend');
pooledf=pooledf(sortedpositions);
osensitivity=cumsum(pooledf>0)./cumsum(ones(size(isolatesf)));
% 
% disp(isolatesf)
%disp(pooledf)

%freqcutoff=find(pooledf <= params.minorfreqthreshold, 1)-1;


% 
% %plot
% figure(88); clf; hold on;
% plot(smooth(ospecificity,10),'r')
% plot(especificity, 'k')
% xlabel('Rank polymoprhic position found in pooled P1')
% ylabel('% of positions polymorphic in deep found in isolates')
% plot([freqcutoff freqcutoff], [0 1]) %bars to indicate cutoff
% 

if figuren > 0 
    
    %plot in terms of absolute
    figure(figuren+1); clf; hold on;
    %background
    freqcutoffs=[12 4 2 1]/29; lowerx=1;
    for i=1:numel(freqcutoffs)
       upperx=find(isolatesf >= freqcutoffs(i), 1, 'last');
       disp([freqcutoffs(i)*29, upperx])
       p1=patch([lowerx; upperx; upperx; lowerx], [0.01; 0.01; numel(isolatesf); numel(isolatesf)],  repmat(1-2*freqcutoffs(i),1,3));
       lowerx=upperx; set(p1,'EdgeColor','None','LineWidth',.01)
    end
    plot(1:numel(isolatesf), 1:numel(isolatesf), 'k:','LineWidth',1.1)
  %  plot(cumsum(isolatesf), 'k', 'LineWidth',4)
    plot(smooth(osensitivity,3),'m', 'LineWidth',4)
    xlabel('Rank polymorphic position in isolates from P1 (log scale)')
    ylabel('Percent of mutations also found pooled sequencing from P1')
   % plot([freqcutoff freqcutoff], [0 max(cumsum(prob_in_isolates))], 'r--','LineWidth',3) %bars to indicate cutoff
    axis([1 upperx 0 1.01]);
    set(gca,'XScale','log')
    set(gca,'Xtick',[1 4 8 16 32 64 128])
% 
%     %plot sensitivity
%     figure(100+ figuren); clf; hold on;
%     freqcutoffs=[.15 .08 .05 .03 .02 .015]; lowerx=0;
%     for i=1:numel(freqcutoffs)
%        upperx=find(pooledf <= freqcutoffs(i), 1)-1;
%        p1=patch([lowerx; upperx; upperx; lowerx], [0.001; 0.001; 2; 2],  repmat(1-5*freqcutoffs(i),1,3));
%        lowerx=upperx; set(p1,'EdgeColor','None','LineWidth',.01)
%     end
%     plot(smooth(ospecificity./especificity,5),'m', 'LineWidth',4)
%     xlabel('Rank polymoprhic position found in pooled P1')
%     ylabel('Specificity observed / Specificiy expected')
%     plot([freqcutoff freqcutoff], [0 1.5],'r--','LineWidth',3) %bars to indicate cutoff
%     axis([0 upperx 0 1.02]);

end


%disp(['O/E, O, E'])
%disp([x(freqcutoff) ospecificity(freqcutoff) especificity(freqcutoff)] )
%disp(sensitivity)

sensitivity=osensitivity(end); 
