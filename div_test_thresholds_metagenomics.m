function [allthresholdsmet] = div_test_thresholds_metagenomics(cnts, p)

%Tami Lieberman April 2012, could be written more matlaby, without going
%one sample at a time


NSamples=size(cnts,3);
NPositions=size(cnts,2);

allthresholdsmet=zeros(NPositions,NSamples);


%assess others, using controlMAF
for sample=1:NSamples;
    v=div_single_sample_test_thresholds(cnts(:,:,sample), p);
    allthresholdsmet(:,sample)=v;
end

 
