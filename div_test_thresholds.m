function [allthresholdsmet] = div_test_thresholds(cnts, p, coveragethresholds)

%Tami Lieberman April 2012, could be written more matlaby, without going
%one sample at a time --- this hooks up nicely into a single sample version
%used during initial data curation


NSamples=size(cnts,3);
NPositions=size(cnts,2);

allthresholdsmet=zeros(NPositions,NSamples);

%do control first
[controlFreq, ~, ~] = div_major_allele_freq(cnts(:,:,1));


%assess others, using controlMAF
for sample=1:NSamples;
    v=div_single_sample_test_thresholds(cnts(:,:,sample), p, controlFreq, coveragethresholds(:,sample));
    allthresholdsmet(:,sample)=v;
end

 
