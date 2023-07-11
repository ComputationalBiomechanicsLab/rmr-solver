function [MAE_sum_CMC, MAE_sum_RMR_GH, MAE_sum_RMR_noGH, CMC_corr, corr_RMR_GH, corr_RMR_noGH]=EMG_CMC_RMR_plot(CMC01_All,CMC02_All,CMC03_All,EMG01_All,EMG02_All,EMG03_All,inds01,indend01,indmax01,inds02,indend02,indmax02,inds03,indend03,indmax03,muscle_CMC,muscleEMG,MVC,muscle_Name, path_RMR_GH, path_RMR_noGH, muscle_RMR, task)

[m1,~]=size(CMC01_All);
[CMC01,CMCpercMot01]=dataNorm(CMC01_All,inds01,indend01,indmax01,muscle_CMC);

[m2,~]=size(CMC02_All);
[CMC02,CMCpercMot02]=dataNorm(CMC02_All,inds02,indend02,indmax02,muscle_CMC);

[m3,~]=size(CMC03_All);
[CMC03,CMCpercMot03]=dataNorm(CMC03_All,inds03,indend03,indmax03,muscle_CMC);

sampleSize_CMC=min([length(CMC01),length(CMC02),length(CMC03)]);

[CMCpercMot01, index_CMC01] = unique(CMCpercMot01);
[CMCpercMot02, index_CMC02] = unique(CMCpercMot02);
[CMCpercMot03, index_CMC03] = unique(CMCpercMot03);

if sampleSize_CMC == length(CMC01)
    target_CMC=CMCpercMot01;
elseif sampleSize_CMC == length(CMC02)
    target_CMC=CMCpercMot02;
else
    target_CMC=CMCpercMot03;
end

CMC01_aln=interp1(CMCpercMot01,CMC01(index_CMC01),target_CMC,'linear');
CMC02_aln=interp1(CMCpercMot02,CMC02(index_CMC02),target_CMC,'linear');
CMC03_aln=interp1(CMCpercMot03,CMC03(index_CMC03),target_CMC,'linear');

CMC01_aln = movmean(CMC01_aln,30);
CMC02_aln = movmean(CMC02_aln,30);
CMC03_aln = movmean(CMC03_aln,30);

[CMCaverage,CMCstdup,CMCstdlow] = meanstd(CMC01_aln,CMC02_aln,CMC03_aln);

fc =100;
fs=1000;
[b_h,a_h]=butter(4,fc/(fs/2),'high');

fs = 1000;
fnyq= fs/2;
fco=2;
[b_l,a_l]=butter(2,fco*1.25/fnyq);

EMG01_raw = EMG01_All(:,muscleEMG)/MVC;
EMG01_filter = EMGfilter(EMG01_raw,a_h,b_h,a_l,b_l,muscle_CMC);
EMG01_filter = EMG01_filter(CMC01_All(1,1)*1000:CMC01_All(m1,1)*1000);
EMGinds01=round(inds01/m1*length(EMG01_filter)); 
EMGindend01=round(indend01/m1*length(EMG01_filter));
EMGindmax01=round(indmax01/m1*length(EMG01_filter));
[EMG01,EMGpercMot01]=EMGnorm(EMG01_filter,EMGinds01,EMGindend01,EMGindmax01);

EMG02_raw = EMG02_All(:,muscleEMG)/MVC;
EMG02_filter = EMGfilter(EMG02_raw,a_h,b_h,a_l,b_l,muscle_CMC);
EMG02_filter = EMG02_filter(CMC02_All(1,1)*1000:CMC02_All(m2,1)*1000);
EMGinds02=round(inds02/m2*length(EMG02_filter)); 
EMGindend02=round(indend02/m2*length(EMG02_filter));
EMGindmax02=round(indmax02/m2*length(EMG02_filter));
[EMG02,EMGpercMot02]=EMGnorm(EMG02_filter,EMGinds02,EMGindend02,EMGindmax02);

EMG03_raw = EMG03_All(:,muscleEMG)/MVC;
EMG03_filter = EMGfilter(EMG03_raw,a_h,b_h,a_l,b_l,muscle_CMC);
EMG03_filter = EMG03_filter(CMC03_All(1,1)*1000:CMC03_All(m3,1)*1000);
EMGinds03=round(inds03/m3*length(EMG03_filter)); 
EMGindend03=round(indend03/m3*length(EMG03_filter));
EMGindmax03=round(indmax03/m3*length(EMG03_filter));
[EMG03,EMGpercMot03]=EMGnorm(EMG03_filter,EMGinds03,EMGindend03,EMGindmax03);

sampleSize_CMC=min([length(EMG01),length(EMG02),length(EMG03)]);
[EMGpercMot01, index_EMG01] = unique(EMGpercMot01);
[EMGpercMot02, index_EMG02] = unique(EMGpercMot02);
[EMGpercMot03, index_EMG03] = unique(EMGpercMot03);

if sampleSize_CMC == length(EMG01)
    target_EMG=EMGpercMot01;
elseif sampleSize_CMC == length(EMG02)
    target_EMG=EMGpercMot02;
else
    target_EMG=EMGpercMot03;
end

EMG01_aln=interp1(EMGpercMot01,EMG01(index_EMG01),target_EMG,'linear');
EMG02_aln=interp1(EMGpercMot02,EMG02(index_EMG02),target_EMG,'linear');
EMG03_aln=interp1(EMGpercMot03,EMG03(index_EMG03),target_EMG,'linear');

EMG01_aln = movmean(EMG01_aln,30);
EMG02_aln = movmean(EMG02_aln,30);
EMG03_aln = movmean(EMG03_aln,30);

[EMGaverage,EMGstdup,EMGstdlow] = meanstd(EMG01_aln,EMG02_aln,EMG03_aln);

if path_RMR_GH
    [RMR_GH_average,RMR_GH_stdup,RMR_GH_stdlow, target_RMR_GH] = evaluate_mean_std_RMR(path_RMR_GH, muscle_RMR, task, 0);
end
if path_RMR_noGH
    [RMR_noGH_average,RMR_noGH_stdup,RMR_noGH_stdlow, target_RMR_noGH] = evaluate_mean_std_RMR(path_RMR_noGH, muscle_RMR, task, 0);
end

fig = figure;
hold on
set(gca,'fontsize',15,'LineWidth',2)

%plotting EMG data
betweenx_emg = [target_EMG,fliplr(target_EMG)];
betweeny_emg = [EMGstdup, fliplr(EMGstdlow)];
h=fill(betweenx_emg,betweeny_emg,[0.8 0.8 0.8]);
set(h,'EdgeColor','none')
plot(target_EMG, EMGstdup,'Color',[0.8 0.8 0.8],'LineWidth',1.5)

% plotting CMC results
betweenx = [target_CMC,fliplr(target_CMC)];
betweeny = [CMCstdup, fliplr(CMCstdlow)];
h=fill(betweenx,betweeny,[1 0.8 1]);
set(h,'EdgeColor','none')
plot(target_CMC, CMCaverage,'m','LineWidth',1.5)

% plotting RMR results (with GH)
if path_RMR_GH
    betweenx = [target_RMR_GH,fliplr(target_RMR_GH)];
    betweeny = [RMR_GH_stdup, fliplr(RMR_GH_stdlow)];
    h=fill(betweenx,betweeny,[0.8 1 0.8]);
    set(h,'EdgeColor','none')
    plot(target_RMR_GH, RMR_GH_average,'g','LineWidth',1.5)
end

% plotting RMR results (without GH)
if path_RMR_noGH
    betweenx = [target_RMR_noGH,fliplr(target_RMR_noGH)];
    betweeny = [RMR_noGH_stdup, fliplr(RMR_noGH_stdlow)];
    h=fill(betweenx,betweeny,[0.8 0.8 1]);
    set(h,'EdgeColor','none')
    plot(target_RMR_noGH, RMR_noGH_average,'Color', [0.3 0.3 1], 'LineWidth',1.5)
end

axis([0 200 0 0.6])
title(muscle_Name)
% optionally save the figures
% figure_name = append(task, '_', muscle_Name{1});
% saveas(fig, figure_name)

% calculate MAE for CMC
sampleSize_CMC=min([length(EMGaverage),length(CMCaverage)]);
if length(EMGaverage)>length(CMCaverage)
    EMG_MAE_1=interp1(target_EMG,EMGaverage,target_CMC,'linear');
    CMC_MAE=CMCaverage;
else
    CMC_MAE=interp1(target_CMC,CMCaverage,target_EMG,'linear');
    EMG_MAE_1=EMGaverage;
end
MAE_diff_CMC=abs(CMC_MAE-EMG_MAE_1);
MAE_sum_CMC=sum(MAE_diff_CMC)/sampleSize_CMC;

if path_RMR_GH
    % calculate MAE for RMR with GH(always less points than the EMG data)
    sampleSize_RMR=min([length(EMGaverage),length(RMR_GH_average)]);
    EMG_MAE_2=interp1(target_EMG,EMGaverage,target_RMR_GH,'linear');
    
    MAE_diff_RMR_GH=abs(RMR_GH_average-EMG_MAE_2);
    MAE_sum_RMR_GH=sum(MAE_diff_RMR_GH)/sampleSize_RMR;
else
    MAE_sum_RMR_GH = "not evaluated";
end

if path_RMR_noGH
    % calculate MAE for RMR without GH(always less points than the EMG data)
    sampleSize_RMR=min([length(EMGaverage),length(RMR_noGH_average)]);
    EMG_MAE_2=interp1(target_EMG,EMGaverage,target_RMR_noGH,'linear');
    
    MAE_diff_RMR_noGH=abs(RMR_noGH_average-EMG_MAE_2);
    MAE_sum_RMR_noGH=sum(MAE_diff_RMR_noGH)/sampleSize_RMR;
else
    MAE_sum_RMR_noGH = "not evaluated";
end

% calculate the cross-correlation for CMC
CMC_corr = xcorr(EMG_MAE_1, CMC_MAE, 0,'coeff');     % retain correlation corresponding to 0 lag, as signals are already aligned

if path_RMR_GH
    % in our case, the minimum dimension will always correspond to the
    % RMR results. Only that case is considered, code can be extended if
    % needed
    EMG_corr=interp1(target_EMG,EMGaverage,target_RMR_GH,'linear');
    corr_RMR_GH = xcorr(EMG_corr, RMR_GH_average, 0, 'coeff');     % retain correlation corresponding to 0 lag, as signals are already aligned
else
    corr_RMR_GH = "not evaluated";
end

if path_RMR_noGH
    % in our case, the minimum dimension will always correspond to the
    % RMR results. Only that case is considered, code can be extended if
    % needed
    EMG_corr=interp1(target_EMG,EMGaverage,target_RMR_noGH,'linear');
    corr_RMR_noGH = xcorr(EMG_corr, RMR_noGH_average, 0, 'coeff');     % retain correlation corresponding to 0 lag, as signals are already aligned
else
    corr_RMR_noGH = "not evaluated";
end