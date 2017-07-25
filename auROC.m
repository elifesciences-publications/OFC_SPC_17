function [AUC] = auROC(basebins,testbins)

% gives area-under-the-curve of a ROC graph
% for the comparison of baseline (basebins)
% spiking to all time periods (testbins) assuming
% binned spiking data (e.g. 100 ms bins)
% as per figure S1 of cohen and uchida 2012 (nature)


%"We compared the histogram of spike counts during the baseline period
% to that during a given bin (solid line) 
%by moving a criterion from zero to the maximum firing rate. 
% We then plotted the probability that the activity during r1 was greater than the
% criteria against the probability that the baseline activity was greater than the criteria. The area under this
% curve quantifies the degree of overlap between the two spike count distributions (i.e., the discriminability
% of the two)."

%% linearize baseline and test bins to single vectors in case they were matrices
basebins=basebins(:);
testbins=testbins(:);

%% normalize all bins (baseline and test) to span 0-to-1 range where 0 == min and 1 == max 
nbasebins = (basebins-min([basebins;testbins]))/(-min([basebins;testbins])+max([basebins;testbins]));
ntestbins = (testbins-min([basebins;testbins]))/(-min([basebins;testbins])+max([basebins;testbins])); 

%% initialize variables for looping through testbin-vs-basebin comparison
count = 0;
basestat=zeros(100,1);
teststat=zeros(100,1);

%% in small steps, calculate x and y axes for ROC plot, though slower with smaller steps

% i.e. for each step, what fraction of baseline bins are greater than a test statistic (ii)? (x-axis)
% these are comparable to 'false positives' in a traditional ROC

% and for each step, what fraction of 'test' bins are greater than a test statistic (ii)? (y-axis)
% these are comparable to 'true positives' in a traditional ROC

for ii = 0:.01:1.01
   count = count+1;
    basestat(count) = sum(nbasebins>=ii)/length(basebins); 
    teststat(count) = sum(ntestbins>=ii)/length(testbins);
    
end

%% plot the ROC curve, if you'd like
% figure;
% subplot(2,1,1)
% hist([nbasebins,ntestbins])
%  subplot(2,1,2)
%  plot(basestat,teststat)

%% calculate the AU of the ROC via simple integration, unless there's missing data
if sum([basebins;testbins])==0
    AUC = NaN;  
else  
AUC = -trapz(basestat,teststat);
end