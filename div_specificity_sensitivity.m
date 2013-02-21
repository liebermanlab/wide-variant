function [specificity, sensitivity, normalized_specificity, expected_sensitivity] = div_specificity_sensitivity(iAF,pAF, controlAF, samplecnts, params, covthresholdssample, majornt, minornt)

%December 2012 TDL
%Modified from div_analyze_sensitivity_using_isolates
%Takes a lot of inputs so as not to repeat many calculations


testpos=div_single_sample_test_thresholds(samplecnts, params, controlAF', covthresholdssample, majornt, minornt);
pAF(testpos==0)=0;

%specificity
isolatesf=iAF(testpos>0);
prob_in_isolates=1.-((1-pAF(testpos>0)).^29);
normalized_specificity=sum(isolatesf>0)/sum(prob_in_isolates);

bothmethods=sum(pAF(iAF>0 & iAF <1)>0);

specificity=bothmethods/sum(testpos);
sensitivity=bothmethods/sum(iAF>0 & iAF <1);

expected_sensitivity=sum(prob_in_isolates)/numel(isolatesf);

