function Go_nogo_SpkTrain_Matrix_generation(inFile,AnalyzeSes,timeBin)

%  Li Yuan, Jan-05-2025
% -------------------------------------------------------------------------
% set parameters for analysis

p.writeToFile = 1;
p.saveDir = 'N:\yuan\AnalysisResults'; % change it to your own folder

p.timeBin = timeBin; % unit sec
p.gauss = p.timeBin*2; % unit sec

% whether read in spike time from whole session or subsessions
% Read in input information
sessInfo = SessInfoImport_Npx(inFile);

for sessInd = AnalyzeSes(1:end)
    close all
    p.delay_len = sessInfo(sessInd).delay_len;
    p.response_len = sessInfo(sessInd).response_len;

    savedir = sprintf('%s%s%s%s%s',p.saveDir,'\ProcessedData\',sessInfo(sessInd).animalID,'-date-',sessInfo(sessInd).recDate);
    if ~exist(savedir, 'dir')
        mkdir(savedir);
    end
    savedir2 = sprintf('%s%s%s%s%s%s',p.saveDir,'\ProcessedData\',sessInfo(sessInd).animalID,'-date-',sessInfo(sessInd).recDate,'\SpkMat_NPY');
    if ~exist(savedir2, 'dir')
        mkdir(savedir2);
    end

    % load corrected TIMESTAMPS
    syncFile = fullfile(savedir,'syncTime.mat');
    syncTime = load(syncFile);
    IO_ts = syncTime.Npx_timeStamps.IO.Ts_Sync;
    clear syncTime

    spkTrainMatrix.animalID = sessInfo(sessInd).animalID;
    spkTrainMatrix.recDate = sessInfo(sessInd).recDate;
    spkTrainMatrix.odorPair1 = sessInfo(sessInd).pair1;
    spkTrainMatrix.odorPair2 = sessInfo(sessInd).pair2;
    spkTrainMatrix.timeBin = p.timeBin;
    spkTrainMatrix.gauss = p.gauss;

    spk_file = sprintf('%s%s%s',savedir,'\npx_Cluster_','imec0');
    load(spk_file); % output:npx_Cluster
    clusterNum = npx_Cluster.clusterNum;

    % read in behavior time
    behaveFile = fullfile(savedir,'behaviorLabel.mat');
    load(behaveFile); % output: behaviorLabel

    if ~isempty(sessInfo(sessInd).pair1)
        trialNum = behaviorLabel.trialInfo_pair1.num;
        onsetInd = behaviorLabel.trialInfo_pair1.onsetInd;
        offsetInd = behaviorLabel.trialInfo_pair1.offsetInd;
        onsetTs = IO_ts(onsetInd);
        offsetTs = IO_ts(offsetInd);
        odor = behaviorLabel.trialInfo_pair1.odor; % 1,2,3,4
        odorNum = length(unique(odor));
        successLabel = behaviorLabel.trialInfo_pair1.success; % 0, premature lick, 1, successful trial
        outcomeLabel = behaviorLabel.trialInfo_pair1.outcome; % 1: Hit, 2: Miss, 3: FalseAlarm, 4: CorrectRejection

        blockStartTime = min(onsetTs);
        blockEndTime = max(offsetTs) + p.delay_len + p.response_len;

        timeBin = blockStartTime:p.timeBin:blockEndTime;
        binNum = length(timeBin) - 1;
        binTrial_IndLabel = zeros(1,binNum);
        binTrial_OdorLabel = zeros(1,binNum);
        binTrial_SuccessLabel = zeros(1,binNum);
        binTrial_OutcomeLabel = zeros(1,binNum);
        spkTrain_Raw = zeros(clusterNum,binNum);
        spkTrain_Gauss = zeros(clusterNum,binNum);

        for n = 1:binNum
            binStart = timeBin(n);
            binEnd = timeBin(n+1);
            idx = arrayfun(@(x) find(onsetTs <= x, 1, 'last'), binStart); % get the nearest past onsetTs
            binTrial_IndLabel(n) = idx;
            binTrial_OdorLabel(n) = odor(idx);
            binTrial_SuccessLabel(n) = successLabel(idx);
            binTrial_OutcomeLabel(n) = outcomeLabel(idx);
        end
        if p.writeToFile
            fileName = sprintf('%s','pair1_Timebin.npy');
            writeNPY(timeBin(1:end-1), fullfile(savedir2, fileName));
            fileName = sprintf('%s','pair1_binTrial_IndLabel.npy');
            writeNPY(binTrial_IndLabel, fullfile(savedir2, fileName));
            fileName = sprintf('%s','pair1_binTrial_OdorLabel.npy');
            writeNPY(binTrial_OdorLabel, fullfile(savedir2, fileName));
            fileName = sprintf('%s','pair1_binTrial_SuccessLabel.npy');
            writeNPY(binTrial_SuccessLabel, fullfile(savedir2, fileName));
            fileName = sprintf('%s','pair1_binTrial_OutcomeLabel.npy');
            writeNPY(binTrial_OutcomeLabel, fullfile(savedir2, fileName));
        end

        for k = 1:clusterNum
            tSpTemp = npx_Cluster.clusterInfo{k}.spkTs;
            fireTemp_Raw = zeros(1,binNum);
            parfor n = 1:binNum
                fireTemp_Raw(n) = sum(tSpTemp>timeBin(n) & tSpTemp<=timeBin(n+1));
            end
            % spkTrain(k,:) = zscore(fireTemp);
            % smooth spike trains to avoid binning time issue
            fireTemp_Gauss = gaussfilt(1:length(fireTemp_Raw),fireTemp_Raw,p.gauss./p.timeBin);
            spkTrain_Raw(k,:) = fireTemp_Raw;
            spkTrain_Gauss(k,:) = fireTemp_Gauss;

        end
        spkTrainMatrix.odorPair1_timeBin = timeBin(1:end-1);
        spkTrainMatrix.odorPair1_timeBin_IndLabel = binTrial_IndLabel;
        spkTrainMatrix.odorPair1_timeBin_OdorLabel = binTrial_OdorLabel;
        spkTrainMatrix.odorPair1_timeBin_SuccessLabel = binTrial_SuccessLabel;
        spkTrainMatrix.odorPair1_timeBin_OutcomeLabel = binTrial_OutcomeLabel;
        spkTrainMatrix.odorPair1_timeBin_spkTrain_Raw = spkTrain_Raw;
        spkTrainMatrix.odorPair1_timeBin_spkTrain_Gauss = spkTrain_Gauss;
    end

    if p.writeToFile
        fileName = sprintf('%s','pair1_spkTrain_Raw.npy');
        writeNPY(spkTrain_Raw, fullfile(savedir2, fileName));
        fileName = sprintf('%s','pair1_spkTrain_Gauss.npy');
        writeNPY(spkTrain_Gauss, fullfile(savedir2, fileName));
        fileName = sprintf('%s%d%s','spkMat-',p.timeBin*10^3,'ms.mat');
        save(fullfile(savedir,fileName), 'spkTrainMatrix','-v7.3');
    end
    clear spkTrainMatrix
    fprintf('Finished spike matrix analysis for session %d\n',sessInd);
end

end
