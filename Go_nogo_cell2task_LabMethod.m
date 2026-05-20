% this code is to read in spike timestamp and align it to the go-nogo event
% time
% Li YUAN, Tohoku, 2026-Mar
% updated to show whole trial ratemap, including ITI
% Li YUAN, Tohoku, 2026-Apr
function Go_nogo_cell2task_LabMethod(inFile,analyzeSes)
close all
% set(0, 'DefaultFigureVisible', 'off');
p.saveFile = 1;
p.savePlot = 1;
p.saveDir = 'N:\yuan\AnalysisResults'; % change it to your own folder

p.bin = 50/10^3; % unit sec, time bin for calculate rate
premsec = 1; % unit sec, time shown before onset of odor
postmsec = 4;
p.odor_len = 1;
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

    goNogo_CellTask_LabMethod.animalID = sessInfo(sessInd).animalID;
    goNogo_CellTask_LabMethod.recDate = sessInfo(sessInd).recDate;
    goNogo_CellTask_LabMethod.odorPair1 = sessInfo(sessInd).pair1;
    goNogo_CellTask_LabMethod.odorPair2 = sessInfo(sessInd).pair2;
    goNogo_CellTask_LabMethod.binWidth = p.bin;


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
            savedir_plot = sprintf('%s%s%s%s%s%s%s',p.saveDir,'\Figures\',sessInfo(sessInd).animalID,'-date-',sessInfo(sessInd).recDate,'\cell2Task_labOrigin\',probeName);
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

        for pairInd = 1:2 % assuming maximum 2 pairs each day
            pairName = sprintf('pair%d',pairInd);
            trialInfoName = sprintf('trialInfo_pair%d',pairInd);

            if ~isempty(sessInfo(sessInd).(pairName))

                odorName = sessInfo(sessInd).(pairName);
                odorTemp = behaviorLabel.(trialInfoName).odor; % 1,2,3,4
                odorNum = length(unique(odorTemp));
                successLabel = behaviorLabel.(trialInfoName).success == 1; % 0, premature lick, 1, successful trial
                outcomeLabel = behaviorLabel.(trialInfoName).outcome; % 1: Hit, 2: Miss, 3: FalseAlarm, 4: CorrectRejection

                % trialNum = behaviorLabel.(trialInfoName).num;
                odor = odorTemp(successLabel);
                onsetInd = behaviorLabel.(trialInfoName).onsetInd(successLabel);
                offsetInd = behaviorLabel.(trialInfoName).offsetInd(successLabel);
                onset_ts = IO_ts(onsetInd);

                goNogo_CellTask_LabMethod.(probeName).(pairName).odor = odor;
                goNogo_CellTask_LabMethod.(probeName).(pairName).onset_ts = onset_ts;

                for k = 1:clusterNum
                    h = figure;
                    h.Position = [100,100,1600,800];

                    tsp = npx_Cluster.clusterInfo{k}.spkTs;
                    xcoord(k) = npx_Cluster.clusterInfo{k}.channel.xcoords;
                    ycoord(k) = npx_Cluster.clusterInfo{k}.channel.ycoords;
                    cluster_ksID(k) = npx_Cluster.clusterInfo{k}.ID;

                    [ax,posAxisPSTH,raster,PSTH_avg,PSTH_se,h2] = gngRaster4odorV2_Nrp(tsp,odor,onset_ts,premsec,postmsec,odorName);

                    if pairInd == 2

                        for axn = 1:odorNum
                            subplot(ax{axn});
                            [b] =  ylim;
                            ylim([0 b(2)]);

                            line([0 0], [[0 b(2)]], 'Color','k','LineStyle','--');
                            line([p.odor_len p.odor_len]*10^3, [[0 b(2)]], 'Color','k','LineStyle','--');
                            line([p.odor_len + p.delay_len p.odor_len + p.delay_len]*10^3, [[0 b(2)]], 'Color','k','LineStyle','--');
                        end

                        if odorNum == 2
                            linkaxes([ax{1},ax{2}],'y');
                            subplot(2,4,1)
                            title( sprintf('%s', 'Early vs Late CORRECT+ERROR'));
                        elseif odorNum == 4
                            linkaxes([ax{1},ax{2}, ax{3},ax{4}],'y');
                            subplot(2,4,1)
                            title( sprintf('%s', 'Familar vs Novel CORRECT+ERROR'));
                        end

                        subplot(2,4,1)
                        title( sprintf('%s', 'Familar vs Novel CORRECT+ERROR'));

                        subplot(2,4,2)
                        cellname = sprintf('%s%s%s%s%d%s%d',sessInfo(sessInd).animalID,'-Day-',sessInfo(sessInd).recDate,'-Depth-',ycoord(k),'-ID-',cluster_ksID(k));
                        title(cellname);

                        if p.savePlot
                            figName = sprintf('%s%s%s%s%s%s%d%s%d%s%s%s',savedir_plot,'\',sessInfo(sessInd).animalID,'-Day-',sessInfo(sessInd).recDate,...
                                '-Depth-',ycoord(k),'-ID-',cluster_ksID(k),'-',pairName,'-cell2Task');
                            print(figName,'-dpng','-r300');
                        end
                        close(h)
                    end

                    goNogo_CellTask_LabMethod.(probeName).(pairName).raster{k} = raster;
                    goNogo_CellTask_LabMethod.(probeName).(pairName).PSTH_avg{k} = PSTH_avg;
                    goNogo_CellTask_LabMethod.(probeName).(pairName).PSTH_se{k} = PSTH_se;
                    goNogo_CellTask_LabMethod.(probeName).(pairName).posAxisPSTH{k} = posAxisPSTH;
                    goNogo_CellTask_LabMethod.(probeName).(pairName).h2{k} = h2;

                end


                if odorNum == 2
                    data_name.raster={'[1]Go_Familiar, [2]Nogo_Familiar'};
                    data_name.PSTH10trials={'[1]First_10trials, [2]Last_10trials,[3]All'};
                    data_name.PSTH5trials={'[1]First_5trials, [2]Last_5trials,[3]All'};
                elseif odorNum == 4
                    data_name.raster={'[1]Go_Familiar, [2]Nogo_Familiar, [3]Go_Novel, [4]Nogo_Novel'};
                    data_name.PSTH10trials={'[1]First_10trials, [2]Last_10trials,[3]All'};
                    data_name.PSTH5trials={'[1]First_5trials, [2]Last_5trials,[3]All'};

                end

                % lab original code, keep it here as I might want to
                % save it in original way   -- Li YUAN
                %
                % task_time.Nlx_O=Nlx_O;
                % task_time.Nlx_TS=Nlx_TS;
                %
                % data_name.raster={'[1]Go_Familiar, [2]Nogo_Familiar, [3]Go_Novel, [4]Nogo_Novel'};
                % data_name.PSTH10trials={'[1]First_10trials, [2]Last_10trials,[3]All'};
                % data_name.PSTH5trials={'[1]First_5trials, [2]Last_5trials,[3]All'};
                % if blockID == blockID_AB
                %     outf=strcat(F{i}(1:end-2),'_PSTH_AB.mat');    %Save TT.PSTH file: note this is for ABCD session in the first block. TODO:rename with a clearer name
                %     save(outf,'raster','posAxisPSTH','PSTH_avg','PSTH_se','data_name','task_time');
                % elseif  blockID == blockID_ABCD
                %     outf=strcat(F{i}(1:end-2),'_PSTH.mat');    %Save TT.PSTH file: note this is for ABCD session in the second block. TODO:rename with a clearer name
                %     save(outf,'raster','posAxisPSTH','PSTH_avg','PSTH_se','data_name','task_time');
                %
                % end
                goNogo_CellTask_LabMethod.(probeName).(pairName).data_name = data_name;

            end
        end

    end


    if p.saveFile
        fileName = 'goNogo_CellTask_LabMethod.mat';
        save(fullfile(savedir,fileName), 'goNogo_CellTask_LabMethod','-v7.3');
    end
    fprintf('Finish behavioral analysis: %d\n', sessInd);
    clear goNogo_CellTask_LabMethod

end
close all
end

