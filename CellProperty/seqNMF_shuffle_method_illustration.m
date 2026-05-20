% this code is to generate the figure for method illustration
% Li Yuan, Apr-2025
close all
clear all

addpath(genpath('C:\Users\Li\Documents\MATLAB\seqNMF-master'));

inFile = 'DualProbe_Fig8Batch_mec_pfc.xlsx';
% Read in input information
sessInfo = SessInfoImport_Dual(inFile);

% which session and seq to plot as example
i = 16; % 1138-20240405
seq_k = 1;
cell_k = [1132,1277];

% Parameters for the SeqNMF
p.timeBin = 300/10^3;
p.L_time = 5; % unit in second
p.seqThresPct = 99;
p.w_min = 1;

spkMatName = 'DelayFire_AutoCorr_RndSpike.mat';
spkMatShuffleFile = fullfile(sessInfo(i).NIDQ,'Processed',spkMatName);
load(spkMatShuffleFile);
p.gaussSigma_time = DelayFire_AutoCorr_RndSpike.gaussSigma_time;
p.timeBin_spkMat = DelayFire_AutoCorr_RndSpike.timeBin;

% load cell int pyr label in MEC
cellTypeFile = fullfile(sessInfo(i).NIDQ,'Processed', 'mec_clusterType.mat');
load(cellTypeFile);

% get valid cell ind
%     pyrLabel_imec1 = mec_clusterType.labelInd;
pyrLabel_imec1 = mec_clusterType.labelInd;
avgRate_imec1 = mec_clusterType.avgRate;
depth_imec1 = mec_clusterType.clusterDepth;
cellLayer_imec1 = mec_clusterType.layerInd;
depthLayer_imec1 = mec_clusterType.posInd;
clusterID_imec1 = mec_clusterType.clusterID;
clusterNum = mec_clusterType.clusterNum;
pyr = (pyrLabel_imec1 == 1)';

cellInd = size(cell_k);
for k = 1:length(cell_k)
    cellInd(k) = find(clusterID_imec1 == cell_k(k));
end

figure
%% plot cell rnd spike generation
subplot(3,4,1)
maxT = 30;
autocorr_Val = DelayFire_AutoCorr_RndSpike.imec1.on30.autocorr(cellInd(1),:);
lagTime = DelayFire_AutoCorr_RndSpike.imec1.on30.lagTime(cellInd(1),:);
plot(lagTime,autocorr_Val,'k');
xlim([-maxT maxT]);
ylim([0 1])

subplot(3,4,2)
% get shuffled spike
shuffleSpikeAll = DelayFire_AutoCorr_RndSpike.imec1.on30.shuffleSpike;
binTime = DelayFire_AutoCorr_RndSpike.imec1.on30.binTime;
spks_match2 = squeeze(shuffleSpikeAll(cellInd(1),1,:));
% calculate the autocorr
spks_match2_Gaussian = gaussfilt(1:length(spks_match2),spks_match2,p.gaussSigma_time./p.timeBin_spkMat);
% calculate the new autocorr
count1 = spks_match2_Gaussian;

[autocorr_Val_shuffle,lag_shuffle]=xcorr(count1,count1,length(count1),'coeff');
plot(lagTime,autocorr_Val_shuffle,'k')
xlim([-maxT maxT]);
ylim([0 1])

% repeat it for second cell
subplot(3,4,3)
maxT = 30;
autocorr_Val = DelayFire_AutoCorr_RndSpike.imec1.on30.autocorr(cellInd(2),:);
lagTime = DelayFire_AutoCorr_RndSpike.imec1.on30.lagTime(cellInd(2),:);
plot(lagTime,autocorr_Val,'k');
xlim([-maxT maxT]);
ylim([0 1])

subplot(3,4,4)
% get shuffled spike
shuffleSpikeAll = DelayFire_AutoCorr_RndSpike.imec1.on30.shuffleSpike;
binTime = DelayFire_AutoCorr_RndSpike.imec1.on30.binTime;
spks_match2 = squeeze(shuffleSpikeAll(cellInd(2),1,:));
% calculate the autocorr
spks_match2_Gaussian = gaussfilt(1:length(spks_match2),spks_match2,p.gaussSigma_time./p.timeBin_spkMat);
% calculate the new autocorr
count1 = spks_match2_Gaussian;

[autocorr_Val_shuffle,lag_shuffle]=xcorr(count1,count1,length(count1),'coeff');
plot(lagTime,autocorr_Val_shuffle,'k')
xlim([-maxT maxT]);
ylim([0 1])

%% plot raw seqNMF and shuffled seqNMF
% load sequences
fileName = sprintf('%s%d%s%d%s','Fig8CenterSeq-',ceil(p.L_time*10^3),'ms-',p.timeBin*10^3,'ms.mat');
seq_file = fullfile(sessInfo(i).NIDQ,'Processed', fileName);
load(seq_file);
% start taking cells from sequences
seqNum = Fig8_seq_Center_detect.K;
seqWght_pyr = Fig8_seq_Center_detect.imec1.zpyr.W;
seqStrength = Fig8_seq_Center_detect.imec1.zpyr.H_alltime;

% load shuffled sequences
fileName = sprintf('%s%d%s%d%s','Fig8CenterSeq-Shuffle-',ceil(p.L_time*10^3),'ms-',p.timeBin*10^3,'ms.mat');
seq_file = fullfile(sessInfo(i).NIDQ,'Processed', fileName);
load(seq_file);

seqNum_Shuffle = Fig8_seq_Center_detect_Shuffle.K;
shuffleTimes = Fig8_seq_Center_detect_Shuffle.shuffleTimes;
shuffleWeight = zeros(sum(pyr),shuffleTimes*seqNum_Shuffle);
for n = 1:shuffleTimes
    w_shuffle = Fig8_seq_Center_detect_Shuffle.imec1.zpyr.W{n};
    w_shuffle2 = sum(w_shuffle,3);
    shuffleWeight(:,(n-1)*seqNum_Shuffle+(1:seqNum_Shuffle)) = w_shuffle2;
end
% get top bound for each ceil
w_thres = prctile(shuffleWeight, p.seqThresPct, 2);
w_thres2 = w_thres;
w_thres2(w_thres < p.w_min) = p.w_min;
    
[max_factor, L_sort, max_sort, hybrid] = helper.ClusterByFactor(seqWght_pyr(:,seq_k,:),1);
indSort_2 = hybrid(:,3);
seqWght_pyr_sort_2 = seqWght_pyr(indSort_2,:,:);
seqWght_pyr_shuffle_sort = squeeze(Fig8_seq_Center_detect_Shuffle.imec1.zpyr.W{1}(indSort_2,:,:));
w_thres_sort = w_thres2(indSort_2);

selectSeqW = squeeze(seqWght_pyr_sort_2(:,seq_k,:));
selectSeqW_sum = sum(selectSeqW,2);
pyrInd = find(pyrLabel_imec1 == 1);
pyrID = clusterID_imec1(pyrInd);
depth_Pyr = depth_imec1(pyrInd); % pyr cell depth in original order
pyrID_sort = pyrID(indSort_2); % sorted cluster ID
pyrIndSort = pyrInd(indSort_2); % sorted cell ind in all cells
depth_Pyr_sort = depth_Pyr(indSort_2); % sorted pyr depth
pyr_sort_Select = find(selectSeqW_sum > w_thres_sort); % cell ind in sorted way
PyrIdx_sort_select = pyrID_sort(pyr_sort_Select); % cluster ID
pyrInd_sort_select = pyrIndSort(pyr_sort_Select); % cell Ind in all cells
depth_Pyr_sort_select = depth_Pyr_sort(pyr_sort_Select);

cellInd2 = size(1:5);
for k = 1:5
    cellInd2(k) = pyr_sort_Select(k+5);
end

subplot(2,5,6)
WPlot(seqWght_pyr_sort_2(:,seq_k,:),1);
TITLE1 = sprintf('%s%d%s%s%s','Seq: ',seq_k);
TITLE2 = 'Sort by all pyr cells';
title({TITLE1;TITLE2},'Interpreter','None')
yl = ylim;
TITLE1 = sprintf('%s%d%s%s%s%d','Rat-',sessInfo(i).ratID,'-Day',sessInfo(i).Date,'Seq-Delay-Pyr: ',k);
title(TITLE1);
        
subplot(2,5,7)
WPlot(seqWght_pyr_shuffle_sort(:,seq_k,:),1);
TITLE1 = sprintf('%s%d%s%s%s','Seq: ',seq_k);
TITLE2 = 'Sort by all pyr cells';
title({TITLE1;TITLE2},'Interpreter','None')
% ylim(yl)

subplot(2,5,8)
k = 4;
shuffleWeight_sort = shuffleWeight(indSort_2,:);
Violin(shuffleWeight_sort(cellInd2(k),:),1,'ShowData',false);
plot([0.5 1.5],[w_thres_sort(cellInd2(k)) w_thres_sort(cellInd2(k))],'r');
plot(1,selectSeqW_sum(cellInd2(k)),'r*')

subplot(2,5,9)
plot(selectSeqW_sum,length(selectSeqW_sum):-1:1)
hold on
plot(w_thres_sort,length(selectSeqW_sum):-1:1)
ylim(yl)
TITLE = sprintf('Pyr cell n = %d',length(selectSeqW_sum));
title(TITLE)

subplot(2,5,10)
WPlot(seqWght_pyr_sort_2(pyr_sort_Select,seq_k,:),1);
TITLE1 = sprintf('%s%d. %d cells','Seq: ',seq_k,length(pyr_sort_Select));
TITLE2 = 'Sort by activated cells';
title({TITLE1;TITLE2},'Interpreter','None')