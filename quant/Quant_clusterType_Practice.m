% fill in the code
% 
function Quant_clusterType_Practice(inFile,AnalyzeSes)

p.saveFile = 1;
p.saveDir = 'N:\yuan\AnalysisResults'; % change it to your own folder

% ISI threshold
% p.isi_Thres = 2/10^3; % unit:sec
p.isi_vio_Thres = 2; % percent
p.spkWidth_Thres = 0.4*10^-3; % second
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

    %% calculate diferent type cluster numbers
    % isi violation ratio, average rate
    % example of two ways
    % single_Cluster = npx_Cluster_type.imec0.isi_vio <= isi_vio_Thres & npx_Cluster_type.imec0.avgRate >= p.avgRate_Thres;
    % single_Cluster = npx_Cluster_type.imec0.single_label == 1 & npx_Cluster_type.imec0.avgRate >= p.avgRate_Thres;

    
    fprintf('Finished analysis for session %d\n',sessInd);
end

%% plot figures
h = figure;
h.Position = [100,100,900,400];

end