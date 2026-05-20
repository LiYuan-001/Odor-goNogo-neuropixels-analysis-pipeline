% This code is to read in extracted clusters and decide whether this is a
% single cluster for the rest analysis
% I separate it from 'Preprocess_Npx_extractClusters.m' because I might set
% different threshold, but I do not want to go through extract spikes each
% time
% Li YUAN, Tohoku, 2026-Apr-20
function Quant_clusterType(inFile,AnalyzeSes)

p.saveFile = 1;
p.saveDir = 'N:\yuan\AnalysisResults'; % change it to your own folder

% ISI threshold
% p.isi_Thres = 2/10^3; % unit:sec
% p.isi_vio_Thres = 2; % percent
% p.spkWidth_Thres = 0.4*10^-3; % second
% p.isiWindow = 0.05; % second
p.avgRate_Thres = 0.01; % Hz

sessNum = length(AnalyzeSes);
clusterNum_all = zeros(sessNum,1);

wide_single_num = zeros(sessNum,1);
narrow_single_num = zeros(sessNum,1);
wide_single_ratio = zeros(sessNum,1);
narrow_single_ratio = zeros(sessNum,1);

wide_multi_num = zeros(sessNum,1);
narrow_multi_num = zeros(sessNum,1);
wide_multi_ratio = zeros(sessNum,1);
narrow_multi_ratio = zeros(sessNum,1);

% Read in input information
sessInfo = SessInfoImport_Npx(inFile);

for i = 1:sessNum
    sessInd = AnalyzeSes(i);
    close all

    savedir = sprintf('%s%s%s%s%s',p.saveDir,'\ProcessedData\',sessInfo(sessInd).animalID,'-date-',sessInfo(sessInd).recDate);

    % get the calculated parameters
    clusterTypeFile = fullfile(savedir,'npx_Cluster_type.mat');
    load(clusterTypeFile);

    clusterNum_all(i) = numel(npx_Cluster_type.imec0.single_label);

    single_Cluster = npx_Cluster_type.imec0.single_label == 1 & npx_Cluster_type.imec0.avgRate >= p.avgRate_Thres;
    multi_Cluster = npx_Cluster_type.imec0.single_label == 0 & npx_Cluster_type.imec0.avgRate >= p.avgRate_Thres;
    wide_Cluster = npx_Cluster_type.imec0.wide_label == 1 & npx_Cluster_type.imec0.avgRate >= p.avgRate_Thres;
    narrow_Cluster = npx_Cluster_type.imec0.wide_label == 0 & npx_Cluster_type.imec0.avgRate >= p.avgRate_Thres;

    wide_single_num(i)  = sum(single_Cluster & wide_Cluster);
    narrow_single_num(i)  = sum(single_Cluster & narrow_Cluster);
    wide_single_ratio(i)  = wide_single_num(i) / clusterNum_all(i);
    narrow_single_ratio(i)  = narrow_single_num(i) / clusterNum_all(i);

    wide_multi_num(i)  = sum(multi_Cluster & wide_Cluster);
    narrow_multi_num(i)  = sum(multi_Cluster & narrow_Cluster);
    wide_multi_ratio(i)  = wide_multi_num(i) / clusterNum_all(i);
    narrow_multi_ratio(i)  = narrow_multi_num(i) / clusterNum_all(i);

    fprintf('Finished analysis for session %d\n',sessInd);
end

h = figure;
h.Position = [100,100,900,400];
subplot(1,2,1)
plot(wide_single_num,'ro-')
hold on
plot(narrow_single_num,'b^-')
title('Single clusters numbers')
legend({'Wide','Narrow'})

subplot(1,2,2)
plot(wide_multi_num,'ro-')
hold on
plot(narrow_multi_num,'b^-')
title('Multi clusters numbers')
legend({'Wide','Narrow'})

end