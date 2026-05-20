% this code is to read in spike timestamp and align it to the go-nogo event
% time
% Li YUAN, Tohoku, 2026-Mar
% updated to show whole trial ratemap, including ITI
% Li YUAN, Tohoku, 2026-Apr
function Go_nogo_cell2taskRateMap(inFile,analyzeSes)
close all
p.saveFile = 1;
p.savePlot = 1;
p.saveDir = 'N:\yuan\AnalysisResults'; % change it to your own folder

p.bin = 100/10^3; % unit sec, time bin for calculate rate
p.smoothTime = p.bin * 2; % unit sec, time bin around central bin
p.preTime = 2; % unit sec, time shown before onset of odor
% p.postTime = 4; % unit sec, time shown after LED + reward period

% Read in input information
sessInfo = SessInfoImport_Npx(inFile);

for sessInd = analyzeSes(1:end)

    close all

    p.delay_len = sessInfo(sessInd).delay_len;
    p.response_len = sessInfo(sessInd).response_len;

    if p.saveFile
        savedir = sprintf('%s%s%s%s%s',p.saveDir,'\ProcessedData\',sessInfo(sessInd).animalID,'-date-',sessInfo(sessInd).recDate);
        if ~exist(savedir, 'dir')
            mkdir(savedir);
        end
    end

    goNogo_CellTaskRateMap.animalID = sessInfo(sessInd).animalID;
    goNogo_CellTaskRateMap.recDate = sessInfo(sessInd).recDate;
    goNogo_CellTaskRateMap.odorPair1 = sessInfo(sessInd).pair1;
    goNogo_CellTaskRateMap.odorPair2 = sessInfo(sessInd).pair2;
    goNogo_CellTaskRateMap.binWidth = p.bin;
    goNogo_CellTaskRateMap.smoothTime = p.smoothTime;


    % load corrected TIMESTAMPS
    syncFile = fullfile(savedir,'syncTime.mat');
    syncTime = load(syncFile);
    IO_ts = syncTime.Npx_timeStamps.IO.Ts_Sync;

    %% start processing each probe
    for m = 1:length(sessInfo(sessInd).probeID)
        probeName =  sprintf('imec%d',sessInfo(sessInd).probeID(m));
        if p.savePlot
            % directory for plot figures
            % generate a folder for each rat eah day under the current folder
            savedir_plot = sprintf('%s%s%s%s%s%s%s',p.saveDir,'\Figures\',sessInfo(sessInd).animalID,'-date-',sessInfo(sessInd).recDate,'\TaskRateMap\',probeName);
            if ~exist(savedir_plot, 'dir')
                mkdir(savedir_plot);
            end
            % delete(strcat(savedir_plot,'\*'));
        end

        % read in behavior time
        behaveFile = fullfile(savedir,'behaviorLabel.mat');
        load(behaveFile); % output: behaviorLabel
        % read in cluster spike time
        spk_file = sprintf('%s%s%s',savedir,'\npx_Cluster_',probeName);
        load(spk_file); % output:npx_Cluster

        clusterNum = npx_Cluster.clusterNum;
        xcoord = nan(clusterNum,1);
        ycoord = nan(clusterNum,1);
        cluster_ksID = nan(clusterNum,1);

        smoothBinNum = ceil(p.smoothTime./p.bin);

        for pairInd = 1:2
            pairName = sprintf('pair%d',pairInd);
            trialInfoName = sprintf('trialInfo_pair%d',pairInd);

            if ~isempty(sessInfo(sessInd).(pairName))

                odor = behaviorLabel.(trialInfoName).odor; % 1,2,3,4
                odorNum = length(unique(odor));
                successLabel = behaviorLabel.(trialInfoName).success; % 0, premature lick, 1, successful trial
                outcomeLabel = behaviorLabel.(trialInfoName).outcome; % 1: Hit, 2: Miss, 3: FalseAlarm, 4: CorrectRejection

                % trialNum = behaviorLabel.(trialInfoName).num;
                onsetInd = behaviorLabel.(trialInfoName).onsetInd;
                offsetInd = behaviorLabel.(trialInfoName).offsetInd;
                minITI = round(min(IO_ts(onsetInd(2:end)) - IO_ts(onsetInd(1:end-1))));
                binNum = ceil((minITI)/p.bin);
                preBinNum = ceil(p.preTime/p.bin);
                timeBin = -p.preTime:p.bin:minITI; % 0 is odor onset
                timeBinNum = length(timeBin)-1;
                if (binNum + preBinNum) ~= timeBinNum
                    error('Bin num mismatch')
                end

                % we can organize into 3d matrix as
                % clusterNum x trial x timebin for one odor. But for easy to
                % use in future, we save each cluster as a 2d matrix as
                % trial x timebin
                
                rateMap_all = cell(clusterNum,odorNum); % all only for trials without pre lickrateMap_correct = cell(clusterNum,odorNum);
                rateMap_correct = cell(clusterNum,odorNum);
                rateMap_error = cell(clusterNum,odorNum);
                rateMap_fail = cell(clusterNum,odorNum);
                spkCount_all = cell(clusterNum,odorNum);
                spkCount_correct = cell(clusterNum,odorNum);
                spkCount_error = cell(clusterNum,odorNum);
                spkCount_fail = cell(clusterNum,odorNum);

                rateMap_all_pre = cell(clusterNum,odorNum);
                rateMap_correct_pre = cell(clusterNum,odorNum);
                rateMap_error_pre = cell(clusterNum,odorNum);
                rateMap_fail_pre = cell(clusterNum,odorNum);
                spkCount_all_pre = cell(clusterNum,odorNum);
                spkCount_correct_pre = cell(clusterNum,odorNum);
                spkCount_error_pre = cell(clusterNum,odorNum);
                spkCount_fail_pre = cell(clusterNum,odorNum);

                for k = 1:clusterNum

                    tsp = npx_Cluster.clusterInfo{k}.spkTs;
                    xcoord(k) = npx_Cluster.clusterInfo{k}.channel.xcoords;
                    ycoord(k) = npx_Cluster.clusterInfo{k}.channel.ycoords;
                    cluster_ksID(k) = npx_Cluster.clusterInfo{k}.ID;

                    for thisOdor = 1:odorNum
                        % successful trials
                        % plot and construct data for odor 1,2,3,4
                        h = figure;
                        h.Position = [100,100,1600,1200];

                        trialNum_thisodor = sum(odor==thisOdor);
                        tiralPlotLim = trialNum_thisodor+3;

                        trialInd_success = (successLabel==1 & odor==thisOdor);
                        trialNum_success = sum(trialInd_success);
                        onsetInd_success = onsetInd(trialInd_success);
                        offsetInd_success = offsetInd(trialInd_success);
                        onsetTs_success = IO_ts(onsetInd_success);
                        offsetTs_success = IO_ts(offsetInd_success);

                        trialInd_correct = (successLabel==1 & odor==thisOdor & (outcomeLabel==1 | outcomeLabel==4));
                        trialNum_correct = sum(trialInd_correct);
                        onsetInd_correct = onsetInd(trialInd_correct);
                        offsetInd_correct = offsetInd(trialInd_correct);
                        onsetTs_correct = IO_ts(onsetInd_correct);
                        offsetTs_correct = IO_ts(offsetInd_correct);

                        trialInd_error = (successLabel==1 & odor==thisOdor & (outcomeLabel==2 | outcomeLabel==3));
                        trialNum_error = sum(trialInd_error);
                        onsetInd_error = onsetInd(trialInd_error);
                        offsetInd_error = offsetInd(trialInd_error);
                        onsetTs_error = IO_ts(onsetInd_error);
                        offsetTs_error = IO_ts(offsetInd_error);

                        trialInd_fail = (successLabel==0 & odor==thisOdor);
                        trialNum_fail = sum(trialInd_fail);
                        onsetInd_fail = onsetInd(trialInd_fail);
                        offsetInd_fail = offsetInd(trialInd_fail);
                        onsetTs_fail = IO_ts(onsetInd_fail);
                        offsetTs_fail = IO_ts(offsetInd_fail);


                        goNogo_CellTaskRateMap.(probeName).(pairName).onsetInd_correct = onsetInd_correct;
                        goNogo_CellTaskRateMap.(probeName).(pairName).offsetInd_correct = offsetInd_correct;
                        goNogo_CellTaskRateMap.(probeName).(pairName).onsetInd_error = onsetInd_error;
                        goNogo_CellTaskRateMap.(probeName).(pairName).offsetInd_error = offsetInd_error;
                        goNogo_CellTaskRateMap.(probeName).(pairName).onsetInd_fail = onsetInd_fail;
                        goNogo_CellTaskRateMap.(probeName).(pairName).offsetInd_fail = offsetInd_fail;

                        count_all = zeros(trialNum_success,timeBinNum); % spike count
                        count_correct = zeros(trialNum_correct,timeBinNum); % spike count
                        count_error = zeros(trialNum_error,timeBinNum); % spike count
                        count_fail = zeros(trialNum_error,timeBinNum); % spike count
                        rateMat_all = zeros(trialNum_success,timeBinNum);
                        rateMat_correct = zeros(trialNum_correct,timeBinNum);
                        rateMat_error = zeros(trialNum_error,timeBinNum);
                        rateMat_fail = zeros(trialNum_fail,timeBinNum);
                        
                        if trialNum_success > 0
                            for n = 1:trialNum_success
                                tsTemp = tsp - onsetTs_success(n);
                                countTemp = histcounts(tsTemp, timeBin);
                                rateTemp = countTemp./p.bin;
                                rateMat_all(n,:) = gaussfilt(1:timeBinNum,rateTemp,smoothBinNum);
                                count_all(n,:) = countTemp;
                            end
                        end


                        % from here I can do it by matching the trial index
                        % I will calculate trial by trial for reading
                        % simplicity
                        % Li Yuan
                        subplot(3,2,[2,4])
                        if trialNum_correct > 0
                            for n = 1:trialNum_correct
                                tsTemp = tsp - onsetTs_correct(n);
                                countTemp = histcounts(tsTemp, timeBin);
                                rateTemp = countTemp./p.bin;
                                rateMat_correct(n,:) = gaussfilt(1:timeBinNum,rateTemp,smoothBinNum);
                                count_correct(n,:) = countTemp;

                                tsTemp2 = tsTemp(tsTemp>=min(timeBin) & tsTemp<=max(timeBin));
                                xPoints = [tsTemp2';tsTemp2'];
                                ypos = n;
                                yPoints = [ypos+zeros(size(tsTemp2'))-0.3;ypos+zeros(size(tsTemp2'))+0.3];

                                if ~isempty(tsTemp2)
                                    plot(xPoints,yPoints,'Color',[0.2,0,0.8]);
                                    hold on
                                end
                            end
                        end

                        if trialNum_error > 0
                            for n = 1:trialNum_error
                                tsTemp = tsp - onsetTs_error(n);
                                countTemp = histcounts(tsTemp, timeBin);
                                rateTemp = countTemp./p.bin;
                                rateMat_error(n,:) = gaussfilt(1:timeBinNum,rateTemp,smoothBinNum);
                                count_error(n,:) = countTemp;
                                tsTemp2 = tsTemp(tsTemp>=min(timeBin) & tsTemp<=max(timeBin));
                                xPoints = [tsTemp2';tsTemp2'];
                                ypos = (n+trialNum_correct+1);
                                yPoints = [ypos+zeros(size(tsTemp2'))-0.3;ypos+zeros(size(tsTemp2'))+0.3];

                                if ~isempty(tsTemp2)
                                    plot(xPoints,yPoints,'R');
                                    hold on
                                end
                            end
                        end

                        if trialNum_fail > 0
                            for n = 1:trialNum_fail
                                tsTemp = tsp - onsetTs_fail(n);
                                countTemp = histcounts(tsTemp, timeBin);
                                rateTemp = countTemp./p.bin;
                                rateMat_fail(n,:) = gaussfilt(1:timeBinNum,rateTemp,smoothBinNum);
                                count_fail(n,:) = countTemp;
                                tsTemp2 = tsTemp(tsTemp>=min(timeBin) & tsTemp<=max(timeBin));

                                xPoints = [tsTemp2';tsTemp2'];
                                ypos = (n+trialNum_error+trialNum_correct+2);
                                yPoints = [ypos+zeros(size(tsTemp2'))-0.3;ypos+zeros(size(tsTemp2'))+0.3];

                                if ~isempty(tsTemp2)
                                    plot(xPoints,yPoints,'k');
                                    hold on
                                end
                            end
                        end


                        % label odor area
                        odorLen = offsetTs_correct(1)-onsetTs_correct(1);
                        patchX = [0,0,odorLen,odorLen];
                        patchY = [0,tiralPlotLim,tiralPlotLim,0];
                        patch(patchX,patchY,'blue','FaceAlpha',.1,'LineStyle','none')
                        % label LED and duration
                        patchX = [odorLen+p.delay_len,odorLen+p.delay_len,odorLen+p.delay_len+p.response_len,odorLen+p.delay_len+p.response_len];
                        patchY = [0,tiralPlotLim,tiralPlotLim,0];
                        patch(patchX,patchY,'y','FaceAlpha',.1,'LineStyle','none')

                        ylim([0 tiralPlotLim])
                        set(gca,'YDir','Reverse')
                        xlim([timeBin(1) timeBin(end)])
                        title('Blue: correct; Red: error; Black: fail')

                        % plot heat map
                        subplot(3,2,[1,3])
                        emptyFill = nan(1,timeBinNum);
                        imageMat = [rateMat_correct;emptyFill;rateMat_error;emptyFill;rateMat_fail];
                        hh=imagesc((1:timeBinNum)*p.bin-p.preTime,1:size(imageMat,1),imageMat);
                        set(hh, 'AlphaData', ~isnan(imageMat));   % NaN becomes transparent
                        set(gca, 'Color', 'w');           % background shows through as white
                        TITLE1 = sprintf('%s%s%s%s%d%s%d',sessInfo(sessInd).animalID,'-Day-',sessInfo(sessInd).recDate,'-Depth-',ycoord(k),'-ID-',cluster_ksID(k));
                        TITLE2 = sprintf('%s-Odor-%s,correct:%d error:%d fail:%d',pairName,sessInfo(sessInd).(pairName){thisOdor},trialNum_correct,trialNum_error,trialNum_fail);
                        title({TITLE1;TITLE2},'Interpreter','None')
                        xlabel('Time')
                        ylabel('Trials')
                        ylim([0 tiralPlotLim])
                        xline(0,'r--')
                        xline(odorLen,'r--')
                        xline(odorLen+p.delay_len,'r--')
                        xline(odorLen+p.delay_len+p.response_len,'r--')
                        xTick= [0 odorLen odorLen+p.delay_len];
                        xTickLabel = {'on','off','LED'};
                        set(gca, 'XTick', xTick, 'XTickLabel', xTickLabel, 'TickLength', [0, 0]);

                        % plot rate matrix
                        subplot(3,2,5)
                        % stdshade(rateMat_fail,0.2,'k',timeBin(1:end-1));
                        hold on
                        if size(rateMat_error,1) >= 5
                            steshade(rateMat_error,0.2,'r',timeBin(1:end-1));
                        end

                        if size(rateMat_correct,1) >= 5
                            steshade(rateMat_correct,0.2,[0.2,0,0.8],timeBin(1:end-1));
                        end
                        xlim([timeBin(1) timeBin(end)])
                        YLIM = ylim;

                        patchX = [0,0,odorLen,odorLen];
                        patchY = [min(YLIM),max(YLIM),max(YLIM),min(YLIM)];
                        patch(patchX,patchY,'blue','FaceAlpha',.1,'LineStyle','none')
                        % label LED and duration
                        patchX = [odorLen+p.delay_len,odorLen+p.delay_len,odorLen+p.delay_len+p.response_len,odorLen+p.delay_len+p.response_len];
                        patchY = [min(YLIM),max(YLIM),max(YLIM),min(YLIM)];
                        patch(patchX,patchY,'y','FaceAlpha',.1,'LineStyle','none')
                        title('Rate with se')
                        ylabel('Hz')

                        figName = sprintf('%s%s%s%s%s%s%d%s%d%s%s%s%s%s',savedir_plot,'\',sessInfo(sessInd).animalID,'-Day-',sessInfo(sessInd).recDate,...
                            '-Depth-',ycoord(k),'-ID-',cluster_ksID(k),'-',pairName,'-',sessInfo(sessInd).(pairName){thisOdor},'-RateMap');
                        print(figName,'-dpng','-r300');
                        close(h)
             
                        rateMap_all{k,thisOdor} = rateMat_all(:,preBinNum+1:end);
                        rateMap_correct{k,thisOdor} = rateMat_correct(:,preBinNum+1:end);
                        rateMap_error{k,thisOdor} = rateMat_error(:,preBinNum+1:end);
                        rateMap_fail{k,thisOdor} = rateMat_fail(:,preBinNum+1:end);
                        spkCount_all{k,thisOdor} = count_all(:,preBinNum+1:end);
                        spkCount_correct{k,thisOdor} = count_correct(:,preBinNum+1:end);
                        spkCount_error{k,thisOdor} = count_error(:,preBinNum+1:end);
                        spkCount_fail{k,thisOdor} = count_fail(:,preBinNum+1:end);

                        rateMap_all_pre{k,thisOdor} = rateMat_all(:,1:preBinNum);
                        rateMap_correct_pre{k,thisOdor} = rateMat_correct(:,1:preBinNum);
                        rateMap_error_pre{k,thisOdor} = rateMat_error(:,1:preBinNum);
                        rateMap_fail_pre{k,thisOdor} = rateMat_fail(:,1:preBinNum);
                        spkCount_all_pre{k,thisOdor} = count_all(:,1:preBinNum);
                        spkCount_correct_pre{k,thisOdor} = count_correct(:,1:preBinNum);
                        spkCount_error_pre{k,thisOdor} = count_error(:,1:preBinNum);
                        spkCount_fail_pre{k,thisOdor} = count_fail(:,1:preBinNum);

                    end
                end

                goNogo_CellTaskRateMap.(probeName).(pairName).rateMap_all = rateMap_all;
                goNogo_CellTaskRateMap.(probeName).(pairName).rateMap_correct = rateMap_correct;
                goNogo_CellTaskRateMap.(probeName).(pairName).rateMap_error = rateMap_error;
                goNogo_CellTaskRateMap.(probeName).(pairName).rateMap_fail = rateMap_fail;
                goNogo_CellTaskRateMap.(probeName).(pairName).spkCount_all = spkCount_all;
                goNogo_CellTaskRateMap.(probeName).(pairName).spkCount_correct = spkCount_correct;
                goNogo_CellTaskRateMap.(probeName).(pairName).spkCount_error = spkCount_error;
                goNogo_CellTaskRateMap.(probeName).(pairName).spkCount_fail = spkCount_fail;

                goNogo_CellTaskRateMap.(probeName).(pairName).rateMap_all_pre = rateMap_all_pre;
                goNogo_CellTaskRateMap.(probeName).(pairName).rateMap_correct_pre = rateMap_correct_pre;
                goNogo_CellTaskRateMap.(probeName).(pairName).rateMap_error_pre = rateMap_error_pre;
                goNogo_CellTaskRateMap.(probeName).(pairName).rateMap_fail_pre = rateMap_fail_pre;
                goNogo_CellTaskRateMap.(probeName).(pairName).spkCount_all_pre = spkCount_all_pre;
                goNogo_CellTaskRateMap.(probeName).(pairName).spkCount_correct_pre = spkCount_correct_pre;
                goNogo_CellTaskRateMap.(probeName).(pairName).spkCount_error_pre = spkCount_error_pre;
                goNogo_CellTaskRateMap.(probeName).(pairName).spkCount_fail_pre = spkCount_fail_pre;
            end
        end
        goNogo_CellTaskRateMap.(probeName).xcoord = xcoord;
        goNogo_CellTaskRateMap.(probeName).ycoord = ycoord;
        goNogo_CellTaskRateMap.(probeName).cluster_ksID = cluster_ksID;
        goNogo_CellTaskRateMap.(probeName).timeBin = timeBin;
        goNogo_CellTaskRateMap.(probeName).preBinNum = preBinNum;
        goNogo_CellTaskRateMap.(probeName).binNum = binNum;
        goNogo_CellTaskRateMap.(probeName).odorLen = odorLen;
        goNogo_CellTaskRateMap.(probeName).delay_len = p.delay_len;
        goNogo_CellTaskRateMap.(probeName).response_len = p.response_len;
    end

    if p.saveFile
        fileName = sprintf('%s%d%s','goNogo_CellTaskRateMap-',ceil(p.bin*1000),'ms.mat');
        save(fullfile(savedir,fileName), 'goNogo_CellTaskRateMap','-v7.3');
    end
    fprintf('Finish behavioral analysis: %d\n', sessInd);
    clear goNogo_CellTaskRateMap

end
close all
end
