function [average,stdup,stdlow] = meanstd(data01,data02,data03)
comp = [data01;data02;data03];
average = mean(comp);stdev = std(comp);
stdup = average+stdev;
stdlow = average-stdev;