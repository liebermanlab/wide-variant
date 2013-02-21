function gcrit = grubbs(values)

%December 2012
%Tami Lieberman 
%Check critical values at http://www.statistics4u.com/fundstat_eng/ee_grubbs_outliertest.html
%gcrit=grubbs(values)

%Performs  grubbs test as described on wikipedia
% http://en.wikipedia.org/wiki/Grubbs'_test_for_outliers

%Returns value of g


gcrit=max(abs(values-mean(values)))/std(values);
disp('Check critical values at http://www.statistics4u.com/fundstat_eng/ee_grubbs_outliertest.html')

