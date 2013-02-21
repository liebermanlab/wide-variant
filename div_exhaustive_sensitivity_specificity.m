function div_exhaustive_sensitivity_specificity


normalized_specificity=div_analyze_specificity_using_isolates(intersect_iMutAF, mutAF, counts, strict_parameters, coveragethresholds, poscontrolsample, 6)
sensitivity=div_analyze_sensitivity_using_isolates(intersect_iMutAF, mutAF, counts, strict_parameters, coveragethresholds, poscontrolsample, 7)


%rocAUC
figure(2);clf;hold on;
plot(0.1:.1:1,0.1:.1:1,':k')
xlabel('1 - Specificity')
ylabel('Sensitivity')

controlaf=.97:.005:.995;
minorfreqthreshold=.015:.005:.035;
bq=17:1:23;
mq=31:1:35;
min_td=5:2:15;
max_td=45:-2:35;
max_sbp=2:4;
max_percent_indels=.05:.05:.2;

k=0;
test_parameters=strict_parameters;
for p1=controlaf
    test_parameters.min_control_MAF=p1;
    for p2=minorfreqthreshold
        test_parameters.minorfreqthreshold=p2;
        for p3=bq
            test_parameters.min_bq=p3;
            for p4=mq
                test_parameters.min_mq=p4;
                for p5=1:numel(min_td)
                    test_parameters.min_td=min_td(p5);
                    test_parameters.max_td=max_td(p5);
                    for p6=max_sbp
                        test_parameters.max_sbp=p6;
                        for p7=max_percent_indels
                            test_parameters.max_percent_indels=p7;

                            [x,y, ~, ~]=div_specificity_sensitivity(intersect_iMutAF, mutAF(:,poscontrolsample), maf(:,1), counts(:,:,poscontrolsample), test_parameters, coveragethresholds(:,poscontrolsample), maNT(:,poscontrolsample)', minorNT(:,poscontrolsample)');
                            % x= 1-  div_single_sample_test_thresholds(cnts(:,:,sample), test_parameters, maf(:,1), coveragethresholds(:,sample));

                            %x=1-div_analyze_specificity_using_isolates(intersect_iMutAF, mutAF, counts, test_parameters, coveragethresholds, poscontrolsample, figuren);

                            %y=div_analyze_sensitivity_using_isolates(intersect_iMutAF, mutAF, counts, test_parameters, coveragethresholds, poscontrolsample, figuren);
                            plot(1-x,y,'.','MarkerSize', 5, 'ButtonDownFcn',{@dispparamters, test_parameters});
                            k=k+1;
                            if floor(k/100)==(k/100); disp(k); end;
                        end
                    end
                end
            end
        end
    end
end

[x,y]=div_specificity_sensitivity(intersect_iMutAF, mutAF(:,poscontrolsample), maf(:,1), counts(:,:,poscontrolsample), strict_parameters, coveragethresholds(:,poscontrolsample), maNT(:,poscontrolsample)', minorNT(:,poscontrolsample)');
plot(1-x,y,'r.','MarkerSize', 5, 'ButtonDownFcn',{@dispparamters, strict_parameters});





%seperate chunk 




controlaf=.97:.005:.995;

figure(4);clf;hold on;
plot(0.1:.1:1,0.1:.1:1,':k')
xlabel('1 -  specificity')
ylabel('Sensitivity')

min_control_MAF=.97:.001:.999;
test_parameters=strict_parameters;
x=[];
y=[];
for i=min_control_MAF
    test_parameters.min_control_MAF=i;
    
    [spec, sens, normalized_spec, expected_spec]=div_specificity_sensitivity(intersect_iMutAF, mutAF(:,poscontrolsample), maf(:,1), counts(:,:,poscontrolsample), test_parameters, coveragethresholds(:,poscontrolsample), maNT(:,poscontrolsample)', minorNT(:,poscontrolsample)');
    plot(1-spec,sens,'k.','MarkerSize', 5, 'ButtonDownFcn',{@dispparamters, test_parameters});
    x(end+1)=1-spec;
    y(end+1)=sens;
end

%plot(x,y,'-')
[specificity, sensitivity, normalized_spec, expected_spec]=div_specificity_sensitivity(intersect_iMutAF, mutAF(:,poscontrolsample), maf(:,1), counts(:,:,poscontrolsample), strict_parameters, coveragethresholds(:,poscontrolsample), maNT(:,poscontrolsample)', minorNT(:,poscontrolsample)');
 plot(1-specificity,sensitivity,'r.','MarkerSize', 5, 'ButtonDownFcn',{@dispparamters, strict_parameters});

