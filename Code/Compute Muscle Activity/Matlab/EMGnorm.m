function [EMGvalue,percMot]=EMGnorm(data,startind,endind,indmax)
up = data(startind:indmax);
perc01 = linspace(0,100,length(up));
down = data(indmax:endind);
percd01 = linspace(100,200,length(down));
percMot = [perc01 percd01];
EMGvalue = [up' down'];