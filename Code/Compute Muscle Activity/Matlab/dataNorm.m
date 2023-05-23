function [dataValue,percMot]=dataNorm(data,startind,endind,indmax,muscle)
up = data(startind:indmax,muscle);
perc01 = linspace(0,100,length(up));
down = data(indmax:endind,muscle);
percd01 = linspace(100,200,length(down));
percMot = [perc01 percd01];
dataValue = [up' down'];