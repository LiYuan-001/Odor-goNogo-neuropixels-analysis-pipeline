% this code is to read in the IO recorded labview data and analyze behavior
% Li YUAN, 2026-Feb
% update to detect manual water trials
% uodate to copy sync file and block file from data folder
% Li YUAN, 2026-May
function Preprocess_Behavior(inFile,analyzeSes)
close all
p.saveFile = 1;
p.savePlot = 1;
p.saveDir = 'N:\yuan\AnalysisResults'; % change it to your own folder

p.dataType = 'A';
p.LC_window = 20; % unit trials
p.LC_window_each = 5;

p.odorThres = 2;
p.LEDThres = 1;
p.lickThres = 1.5;
p.rewardThres = 2;

cTablegreen = [3/255 175/255 122/255];
cTableOrange = [246/255 170/255 0/255];

% Read in input information
sessInfo = SessInfoImport_Npx(inFile);

for sessInd = analyzeSes(1:end)

    close all

    if p.saveFile
        savedir = sprintf('%s%s%s%s%s',p.saveDir,'\ProcessedData\',sessInfo(sessInd).animalID,'-date-',sessInfo(sessInd).recDate); 
        if ~exist(savedir, 'dir')
            mkdir(savedir);
        end
    end

    if p.savePlot
        savedir_plot = sprintf('%s%s%s%s%s%s',p.saveDir,'\Figures\',sessInfo(sessInd).animalID,'-date-',sessInfo(sessInd).recDate,'\LearningCurve');
        if ~exist(savedir_plot, 'dir')
            mkdir(savedir_plot);
        end
    end

    behaviorLabel.animalID = sessInfo(sessInd).animalID;
    behaviorLabel.recDate = sessInfo(sessInd).recDate;
    behaviorLabel.odorPair1 = sessInfo(sessInd).pair1;
    behaviorLabel.odorPair2 = sessInfo(sessInd).pair2;

    % copy the files to result folder if doesnt exist
    syncFileName = 'syncTime.mat';
    syncDataFile   = fullfile(sessInfo(sessInd).npx_path, syncFileName);
    syncResultFile = fullfile(savedir, syncFileName);
    if exist(syncDataFile, 'file') && ~exist(syncResultFile, 'file')
        copyfile(syncDataFile, syncResultFile);
        fprintf('Copied %s to result folder.\n', syncFileName);
    end
    blockFileName = 'behaviorBlocks.mat';
    blockDataFile   = fullfile(sessInfo(sessInd).npx_path, blockFileName);
    blockResultFile = fullfile(savedir, blockFileName);
    if exist(blockDataFile, 'file') && ~exist(blockResultFile, 'file')
        copyfile(blockDataFile, blockResultFile);
        fprintf('Copied %s to result folder.\n', blockFileName);
    end

    % detect if the recording is using onebox or NIDQ
    NIDQMeta = dir(fullfile(sessInfo(sessInd).npx_path, '*NIDQ.meta'));
    obxMeta  = dir(fullfile(sessInfo(sessInd).npx_path, '*obx.meta'));

    % Check if either file exists
    if ~isempty(NIDQMeta)
        fprintf('Found IO.meta file: %s\n', NIDQMeta(1).name);
        IOFile = dir(fullfile(sessInfo(sessInd).npx_path,'*.nidq.bin'));
        metaFile = dir(fullfile(sessInfo(sessInd).npx_path,'*.nidq.meta'));
        IOFileName = IOFile.name;
        metaFileName = metaFile.name;
        IO_type = 'nidq';
    elseif ~isempty(obxMeta)
        fprintf('Found obx.meta file: %s\n', obxMeta(1).name);
        IOFile = dir(fullfile(sessInfo(sessInd).npx_path,'*.obx.bin'));
        metaFile = dir(fullfile(sessInfo(sessInd).npx_path,'*.obx.meta'));
        IOFileName = IOFile.name;
        metaFileName = metaFile.name;
        IO_type = 'obx';
    else
        IO_type = [];
        fprintf('Error,No IO.meta or obx.meta found in %s\n', folderPath);
        return
    end

    %% read event from IO box if file exist
    if ~isempty(IO_type)
        % Parse the corresponding metafile
        meta_IO = ReadMeta(fullfile(sessInfo(sessInd).npx_path,metaFileName));
        IO_Fs = SampRate(meta_IO);
        dataArray = ReadBin(fullfile(sessInfo(sessInd).npx_path,IOFileName),meta_IO);

        if contains(IO_type,'nidq')
            if p.dataType == 'A'
                behaveData = GainCorrectNI(dataArray, [0:7]+1, meta_IO); % change the [0:7] is only part of channels activated in setting
            else
                % (1-based). For imec data there is never more than one saved digital word.
                dw = 1;
                % Read these lines in dw (0-based).
                dLineList = [0:7] + 1; % correct counts for matlab
                digArray = ExtractDigital(dataArray, meta_IO, dw, dLineList);
                behaveData = double(digArray);
                clear digArray
            end

        elseif contains(IO_type,'obx')
            if p.dataType == 'A'
                behaveData = GainCorrectOBX(dataArray, [0:11]+1, meta_IO); % change the [0:11] is only part of channels activated in setting
            else
                % (1-based). For imec data there is never more than one saved digital word.
                dw = 1;
                % Read these lines in dw (0-based).
                dLineList = [0:11] + 1; % correct counts for matlab
                digArray = ExtractDigital(dataArray, meta_IO, dw, dLineList);
                behaveData = double(digArray);
                clear digArray
            end
        end
        clear dataArray

        if exist(syncResultFile, 'file')
            syncTime = load(syncResultFile);
            IO_ts = syncTime.Npx_timeStamps.IO.Ts_Sync;
        else
            IO_ts = (0:size(behaveData,2)-1) / IO_Fs;
        end


        % go_nogo_blockCount = 0;
        % get the split information
        if exist(fullfile(savedir,'behaviorBlocks.mat'),'file')
            load(fullfile(savedir,'behaviorBlocks.mat'));
            eventNames = behaviorBlocks.eventNames;
            eventInd = behaviorBlocks.signalInd;

            if any(ismember(eventNames,'pair1start'))
                go_nogo_blockCount = 1;
                pair1StartInd = eventInd(ismember(eventNames,'pair1start'));
                pair1EndInd = eventInd(ismember(eventNames,'pair1end'));
                behaviorLabel.pair1StartInd = pair1StartInd;
                behaviorLabel.pair1EndInd = pair1EndInd;
            end

            if any(ismember(eventNames,'pair2start'))
                go_nogo_blockCount = 2;
                pair2StartInd = eventInd(ismember(eventNames,'pair2start'));
                pair2EndInd = eventInd(ismember(eventNames,'pair2end'));
                behaviorLabel.pair2StartInd = pair2StartInd;
                behaviorLabel.pair2EndInd = pair2EndInd;
            end

        else
            go_nogo_blockCount = 1;
            pair1StartInd = 1;
            pair1EndInd = size(behaveData,2);
        end

        % get each signal line
        if go_nogo_blockCount == 2

            behaveData_pair1 = behaveData(:,pair1StartInd:pair1EndInd);
            behaveData_pair2 = behaveData(:,pair2StartInd:pair2EndInd);

            % detect the trial onset index, trial related info (odor,
            % outcome etc)
            trialInfo_pair1 = taskPerformance_onebox(behaveData_pair1,IO_Fs,sessInfo(sessInd),'pair1',p);
            trialInfo_pair2 = taskPerformance_onebox(behaveData_pair2,IO_Fs,sessInfo(sessInd),'pair2',p);
            trialInfo_pair1.onsetInd = trialInfo_pair1.onsetInd + pair1StartInd -1;
            trialInfo_pair1.offsetInd = trialInfo_pair1.offsetInd + pair1StartInd -1;
            trialInfo_pair2.onsetInd = trialInfo_pair2.onsetInd + pair2StartInd -1;
            trialInfo_pair2.offsetInd = trialInfo_pair2.offsetInd + pair2StartInd -1;
            trialInfo_pair1.onsetTs = IO_ts(trialInfo_pair1.onsetInd);
            trialInfo_pair1.offsetTs = IO_ts(trialInfo_pair1.offsetInd);
            trialInfo_pair2.onsetTs = IO_ts(trialInfo_pair2.onsetInd);
            trialInfo_pair2.offsetTs = IO_ts(trialInfo_pair2.offsetInd);

            matureTrialAll_pair1 = sum(trialInfo_pair1.mature == 1);
            matureTrialA_pair1 = sum(trialInfo_pair1.mature == 1 & trialInfo_pair1.odor == 1);
            matureTrialB_pair1 = sum(trialInfo_pair1.mature == 1 & trialInfo_pair1.odor == 2);
            matureTrialC_pair1 = sum(trialInfo_pair1.mature == 1 & trialInfo_pair1.odor == 3);
            matureTrialD_pair1 = sum(trialInfo_pair1.mature == 1 & trialInfo_pair1.odor == 4);

            matureTrialAll_pair2 = sum(trialInfo_pair2.mature == 1);
            matureTrialA_pair2 = sum(trialInfo_pair2.mature == 1 & trialInfo_pair2.odor == 1);
            matureTrialB_pair2 = sum(trialInfo_pair2.mature == 1 & trialInfo_pair2.odor == 2);
            matureTrialC_pair2 = sum(trialInfo_pair2.mature == 1 & trialInfo_pair2.odor == 3);
            matureTrialD_pair2 = sum(trialInfo_pair2.mature == 1 & trialInfo_pair2.odor == 4);

            % initiate manual water trial label
            trialInfo_pair1.manual = nan(trialInfo_pair1.num,1);  
            trialInfo_pair1.manual(trialInfo_pair1.mature == 1) = 0;
            trialInfo_pair2.manual = nan(trialInfo_pair2.num,1);
            trialInfo_pair2.manual(trialInfo_pair2.mature == 1) = 0;

            % assign manual and compare trial detection from labview event log
            % compare the two event streams (onebox and labview log)
            % event log can be both longer or shorter than onebox pair
            % length, we sometimes stop labview or cut pair
            % read event from eventLOG.txt file
            log_folder_pair1 = dir(fullfile(sessInfo(sessInd).npx_path, '*pair1'));
            log_folder_pair2 = dir(fullfile(sessInfo(sessInd).npx_path, '*pair2'));
            if ~isempty(log_folder_pair1)
                logFile_pair1 = dir(fullfile(log_folder_pair1.folder,log_folder_pair1.name,'*eventLOG.txt'));
                log_pair1 = load(fullfile(logFile_pair1.folder,logFile_pair1.name));
                [evnt_pair1,results_pair1] = gngGetTS_check_error_AB_ABCD_labView(log_pair1);
                matureInd = trialInfo_pair1.mature == 1;


                % find matching trials
                match = findContinuousOdorMatch(trialInfo_pair1, evnt_pair1);
                matchPct = match.bestLen / sum(matureInd) * 100;
                fprintf('pair1 %2.1f trials found in lab view eventLog\n',matchPct);
                
                oneBoxOutcome = trialInfo_pair1.outcome(matureInd);
                oneBoxOutcome2 = oneBoxOutcome(match.oneboxStart:match.oneboxEnd);
                labViewOutcome = evnt_pair1.action(match.eventStart:match.eventEnd);
                oneBoxOutcome2 = oneBoxOutcome2(:);
                labViewOutcome = labViewOutcome(:);

                oneBoxManual = zeros(sum(matureInd),1);
                oneBoxManual(match.oneboxStart:match.oneboxEnd) = evnt_pair1.MW(match.eventStart:match.eventEnd);
                trialInfo_pair1.manual(matureInd) = oneBoxManual;

                % for manual water trials, always assign outcome as 2
                % or 4, dependeing on odors
                % oneBoxOdor = trialInfo_pair2.odor(matureInd);
                % oneBoxOdor2 = oneBoxOdor(match.oneboxStart:match.oneboxEnd);
                % manualInd = oneBoxManual == 1;
                % odor 1 or 3 + manual water = outcome 2
                % oneBoxOutcome2(manualInd & ismember(oneBoxOdor2, [1 3])) = 2;
                % % odor 2 or 4 + manual water = outcome 3
                % oneBoxOutcome2(manualInd & ismember(oneBoxOdor2, [2 4])) = 3;

                % Assign back
                % oneBoxOutcome(match.oneboxStart:match.oneboxEnd) = oneBoxOutcome2;
                % trialInfo_pair2.outcome(matureInd) = oneBoxOutcome;
                % 
                % if any(oneBoxOutcome2(:) ~= labViewOutcome(:))
                %     error('Pair 2 Outcome mismatch between OneBox and LabView eventLog.');
                % else
                %     fprintf('pair2 %2.1f trials found in lab view eventLog\n',matchPct);
                %     fprintf('pair2 all trial outcome matches between onebox and labview\n');
                % end

            else
                warning('No labview_pair1 event folder')
            end

            if ~isempty(log_folder_pair2)
                logFile_pair2 = dir(fullfile(log_folder_pair2.folder,log_folder_pair2.name,'*eventLOG.txt'));
                log_pair2 = load(fullfile(logFile_pair2.folder,logFile_pair2.name));
                [evnt_pair2,results_pair2] = gngGetTS_check_error_AB_ABCD_labView(log_pair2);
                matureInd = trialInfo_pair2.mature == 1;


                % find matching trials
                match = findContinuousOdorMatch(trialInfo_pair2, evnt_pair2);
                matchPct = match.bestLen / sum(matureInd) * 100;
                fprintf('pair2 %2.1f trials found in lab view eventLog\n',matchPct);

                oneBoxOutcome = trialInfo_pair2.outcome(matureInd);
                oneBoxOutcome2 = oneBoxOutcome(match.oneboxStart:match.oneboxEnd);
                labViewOutcome = evnt_pair2.action(match.eventStart:match.eventEnd);
                oneBoxOutcome2 = oneBoxOutcome2(:);
                labViewOutcome = labViewOutcome(:);

                oneBoxManual = zeros(sum(matureInd),1);
                oneBoxManual(match.oneboxStart:match.oneboxEnd) = evnt_pair2.MW(match.eventStart:match.eventEnd);
                trialInfo_pair2.manual(matureInd) = oneBoxManual;

                % for manual water trials, always assign outcome as 2
                % or 4, dependeing on odors
                % oneBoxOdor = trialInfo_pair2.odor(matureInd);
                % oneBoxOdor2 = oneBoxOdor(match.oneboxStart:match.oneboxEnd);
                % manualInd = oneBoxManual == 1;
                % odor 1 or 3 + manual water = outcome 2
                % oneBoxOutcome2(manualInd & ismember(oneBoxOdor2, [1 3])) = 2;
                % % odor 2 or 4 + manual water = outcome 3
                % oneBoxOutcome2(manualInd & ismember(oneBoxOdor2, [2 4])) = 3;

                % Assign back
                % oneBoxOutcome(match.oneboxStart:match.oneboxEnd) = oneBoxOutcome2;
                % trialInfo_pair2.outcome(matureInd) = oneBoxOutcome;
                % 
                % if any(oneBoxOutcome2(:) ~= labViewOutcome(:))
                %     error('Pair 2 Outcome mismatch between OneBox and LabView eventLog.');
                % else
                %     fprintf('pair2 %2.1f trials found in lab view eventLog\n',matchPct);
                %     fprintf('pair2 all trial outcome matches between onebox and labview\n');
                % end

            else
                warning('No labview_pair2 event folder')
            end
            behaviorLabel.trialInfo_pair1 = trialInfo_pair1;
            behaviorLabel.trialInfo_pair2 = trialInfo_pair2;

            h = figure(1);
            h.Position = [100,100,1200,900];
            subplot(4,1,1)
            f(1) = subplot(4,1,1);
            ts = IO_ts(pair1StartInd:pair1EndInd);
            plot(ts,behaveData_pair1(sessInfo(sessInd).odorA_ch,:))
            hold on
            plot(ts,behaveData_pair1(sessInfo(sessInd).odorB_ch,:))
            plot(ts,behaveData_pair1(sessInfo(sessInd).odorC_ch,:))
            plot(ts,behaveData_pair1(sessInfo(sessInd).odorD_ch,:))
            % plot(ts,data(sessInfo(sessInd).odorE_ch,:))
            % plot(ts,data(sessInfo(sessInd).odorF_ch,:))
            plot([ts(1),ts(end)],[p.odorThres,p.odorThres],'k')
            text(ts(1),p.odorThres+0.5,'Thres')
            legend({'A','B','C','D'})
            ylabel('Odor Amp (V)')
            TITLE1 = sprintf('%s-%s-pair1-onebox-IO-signal',behaviorLabel.animalID,behaviorLabel.recDate);
            title(TITLE1,'Interpreter','None');


            f(2) = subplot(4,1,2);
            plot(ts,behaveData_pair1(sessInfo(sessInd).led_ch,:))
            hold on
            plot([ts(1),ts(end)],[p.LEDThres,p.LEDThres],'k')
            text(ts(1),p.LEDThres+0.5,'Thres')
            ylabel('LED Amp (V)')
            ylim([0 4])

            f(3) = subplot(4,1,3);
            plot(ts,behaveData_pair1(sessInfo(sessInd).lick_ch,:))
            hold on
            plot([ts(1),ts(end)],[p.lickThres,p.lickThres],'k')
            text(ts(1),p.lickThres+0.5,'Thres')
            % label manual trials
            onsetInd = trialInfo_pair1.onsetInd;
            maualInd = onsetInd(trialInfo_pair1.manual == 1);
            if ~isempty(maualInd)
                plot(IO_ts(maualInd), 4.5*ones(size(ts(maualInd))),'r*')
            end
            ylabel('Lick Amp (V)')
            title('Lick signal not used in deciding correct or not, manual water *')
            ylim([0 5])

            f(4) = subplot(4,1,4);
            plot(ts,behaveData_pair1(sessInfo(sessInd).sucrose_ch,:))
            hold on
            plot(ts,behaveData_pair1(sessInfo(sessInd).quinine_ch,:))
            plot([ts(1),ts(end)],[p.rewardThres,p.rewardThres],'k')
            text(ts(1),p.rewardThres+0.5,'Thres')
            legend({'S','Q'})
            ylabel('Reward Amp (V)')
            xlabel('Time (s)')

            linkaxes(f,'x')
            figName = sprintf('%s%s%s-%s%s',savedir_plot,'\',behaviorLabel.animalID,behaviorLabel.recDate,'-pair1-onebox-IO-signal');
            print(figName,'-dpng','-r200');


            h = figure(2);
            h.Position = [100,100,1200,900];
            subplot(4,1,1)
            f(1) = subplot(4,1,1);
            ts = IO_ts(pair2StartInd:pair2EndInd);
            plot(ts,behaveData_pair2(sessInfo(sessInd).odorA_ch,:))
            hold on
            plot(ts,behaveData_pair2(sessInfo(sessInd).odorB_ch,:))
            plot(ts,behaveData_pair2(sessInfo(sessInd).odorC_ch,:))
            plot(ts,behaveData_pair2(sessInfo(sessInd).odorD_ch,:))
            % plot(ts,data(sessInfo(sessInd).odorE_ch,:))
            % plot(ts,data(sessInfo(sessInd).odorF_ch,:))
            plot([ts(1),ts(end)],[p.odorThres,p.odorThres],'k')
            text(ts(1),p.odorThres+0.5,'Thres')
            legend({'A','B','C','D'})
            ylabel('Odor Amp (V)')
            TITLE1 = sprintf('%s-%s-pair2-onebox-IO-signal',behaviorLabel.animalID,behaviorLabel.recDate);
            title(TITLE1,'Interpreter','None');


            f(2) = subplot(4,1,2);
            plot(ts,behaveData_pair2(sessInfo(sessInd).led_ch,:))
            hold on
            plot([ts(1),ts(end)],[p.LEDThres,p.LEDThres],'k')
            text(ts(1),p.LEDThres+0.5,'Thres')
            ylabel('LED Amp (V)')
            ylim([0 4])

            f(3) = subplot(4,1,3);
            plot(ts,behaveData_pair2(sessInfo(sessInd).lick_ch,:))
            hold on
            plot([ts(1),ts(end)],[p.lickThres,p.lickThres],'k')
            text(ts(1),p.lickThres+0.5,'Thres')
            ylabel('Lick Amp (V)')
            title('Lick signal not used in deciding correct or not, manual water *')
            % label manual trials
            onsetInd = trialInfo_pair2.onsetInd;
            maualInd = onsetInd(trialInfo_pair2.manual == 1);
            if ~isempty(maualInd)
                plot(IO_ts(maualInd), 4.5*ones(size(ts(maualInd))),'r*')
            end
            ylim([0 5])

            f(4) = subplot(4,1,4);
            plot(ts,behaveData_pair2(sessInfo(sessInd).sucrose_ch,:))
            hold on
            plot(ts,behaveData_pair2(sessInfo(sessInd).quinine_ch,:))
            plot([ts(1),ts(end)],[p.rewardThres,p.rewardThres],'k')
            text(ts(1),p.rewardThres+0.5,'Thres')
            legend({'S','Q'})
            ylabel('Reward Amp (V)')
            xlabel('Time (s)')
            
            linkaxes(f,'x')
            figName = sprintf('%s%s%s-%s%s',savedir_plot,'\',behaviorLabel.animalID,behaviorLabel.recDate,'-pair2-onebox-IO-signal');
            print(figName,'-dpng','-r200');


            % plot the accuracy plot
            % plot the learning curve for pair 1
            h = figure(3);
            h.Position = [100,100,1200,900];
            subplot(2,1,1)
            [LC,~] = learningCurveCalc(behaviorLabel.trialInfo_pair1,p.LC_window,p.LC_window_each);

            subplot(2,2,1)
            % AB vs CD
            plot(LC.ABCD(1,:),'b','LineWidth',2);hold on;axis tight
            plot(LC.ABCD(4,:),'r','LineWidth',2)
            % plot(matureRate.all,'k','LineWidth',2)
            ylim([30 100])
            label = ['Trials (mov window = ' num2str(p.LC_window) ' trials)'];
            xlabel(label)
            ylabel('Correct %')
            legend('Correct% A/B ','Correct% novel','Location','southeast')
            TITLE1 = sprintf('%s-%s-%s',behaviorLabel.animalID,behaviorLabel.recDate,cell2mat(behaviorLabel.odorPair1));
            TITLE2 = sprintf('Success trial num: %d',matureTrialAll_pair1);
            title({TITLE1;TITLE2},'interpreter','None')


            subplot(2,2,2)
            % A, B, C, D in all type trials
            plot(LC.ABCD(5,:),'Color',cTablegreen,'LineWidth',2);hold on;axis tight
            plot(LC.ABCD(6,:),'Color',cTableOrange,'LineWidth',2)
            plot(LC.ABCD(2,:),'r','LineWidth',2)
            plot(LC.ABCD(3,:),'b','LineWidth',2)
            ylim([0 100])
            label = ['Trials (mov window = ' num2str(p.LC_window) ' trials)'];
            xlabel(label)
            ylabel('Correct %')
            title('Each odor')
            legend('1','2','A','B','Location','southeast')

            subplot(2,2,3)
            plot(LC.C,'Color',cTablegreen,'LineWidth',2);
            hold on
            axis tight
            plot(LC.D,'Color',cTableOrange,'LineWidth',2);
            plot(LC.A,'r','LineWidth',2);
            plot(LC.B,'b','LineWidth',2);
            label = ['Trials (each odor); movWin = ' num2str(p.LC_window_each) ' trials'];
            xlabel(label)
            legend('1','2','A','B','Location','southeast')
            ylim([0 100])
            TITLE1 = sprintf('Trial num A: %d B: %d',matureTrialA_pair1,matureTrialB_pair1);
            TITLE2 = sprintf('Trial num 1: %d 2: %d',matureTrialC_pair1,matureTrialD_pair1);
            title({TITLE1;TITLE2},'interpreter','None')
            
            subplot(2,2,4)
            plot(LC.CD,'r','LineWidth',2);hold on;axis tight
            plot(LC.AB,'b','LineWidth',2);
            label = ['Trials (each odor); movWin = ' num2str(p.LC_window_each) ' trials'];
            xlabel(label)
            legend('12','AB','Location','southeast')
            ylim([0 100])
            TITLE1 = sprintf('Trial num AB %d',matureTrialA_pair1+matureTrialB_pair1);
            TITLE2 = sprintf('Trial num 1-2: %d',matureTrialC_pair1+matureTrialD_pair1);
            title({TITLE1;TITLE2},'interpreter','None')
            figName = sprintf('%s%s%s-%s%s',savedir_plot,'\',behaviorLabel.animalID,behaviorLabel.recDate,'-pair1-learningCurve');
            print(figName,'-dpng','-r200');

            behaviorLabel.LC_pair1 = LC;


            % plot the learning curve for pair 2
            h = figure(4);
            h.Position = [100,100,1200,900];
            subplot(2,1,1)
            [LC,matureRate] = learningCurveCalc(behaviorLabel.trialInfo_pair2,p.LC_window,p.LC_window_each);

            subplot(2,2,1)
            % AB vs CD
            plot(LC.ABCD(1,:),'b','LineWidth',2);hold on;axis tight
            plot(LC.ABCD(4,:),'r','LineWidth',2)
            % plot(matureRate.all,'k','LineWidth',2)
            ylim([30 100])
            label = ['Trials (mov window = ' num2str(p.LC_window) ' trials)'];
            xlabel(label)
            ylabel('Correct %')
            legend('Correct% A/B ','Correct% novel','Location','southeast')
            TITLE1 = sprintf('%s-%s-%s',behaviorLabel.animalID,behaviorLabel.recDate,cell2mat(behaviorLabel.odorPair2));
            TITLE2 = sprintf('Success trial num: %d',matureTrialAll_pair2);
            title({TITLE1;TITLE2},'interpreter','None')


            subplot(2,2,2)
            % A, B, C, D in all type trials
            plot(LC.ABCD(5,:),'Color',cTablegreen,'LineWidth',2);hold on;axis tight
            plot(LC.ABCD(6,:),'Color',cTableOrange,'LineWidth',2)
            plot(LC.ABCD(2,:),'r','LineWidth',2)
            plot(LC.ABCD(3,:),'b','LineWidth',2)
            ylim([0 100])
            label = ['Trials (mov window = ' num2str(p.LC_window) ' trials)'];
            xlabel(label)
            ylabel('Correct %')
            title('Each odor')
            legend('1','2','A','B','Location','southeast')

            subplot(2,2,3)
            plot(LC.C,'Color',cTablegreen,'LineWidth',2);
            hold on
            axis tight
            plot(LC.D,'Color',cTableOrange,'LineWidth',2);
            plot(LC.A,'r','LineWidth',2);
            plot(LC.B,'b','LineWidth',2);
            label = ['Trials (each odor); movWin = ' num2str(p.LC_window_each) ' trials'];
            xlabel(label)
            legend('1','2','A','B','Location','southeast')
            ylim([0 100])
            TITLE1 = sprintf('Trial num A: %d B: %d',matureTrialA_pair2,matureTrialB_pair2);
            TITLE2 = sprintf('Trial num 1: %d 2: %d',matureTrialC_pair2,matureTrialD_pair2);
            title({TITLE1;TITLE2},'interpreter','None')
            
            subplot(2,2,4)
            plot(LC.CD,'r','LineWidth',2);hold on;axis tight
            plot(LC.AB,'b','LineWidth',2);
            label = ['Trials (each odor); movWin = ' num2str(p.LC_window_each) ' trials'];
            xlabel(label)
            legend('12','AB','Location','southeast')
            ylim([0 100])
            TITLE1 = sprintf('Trial num AB %d',matureTrialA_pair2+matureTrialB_pair2);
            TITLE2 = sprintf('Trial num 1-2: %d',matureTrialC_pair2+matureTrialD_pair2);
            title({TITLE1;TITLE2},'interpreter','None')
            figName = sprintf('%s%s%s-%s%s',savedir_plot,'\',behaviorLabel.animalID,behaviorLabel.recDate,'-pair2-learningCurve');
            print(figName,'-dpng','-r200');

            behaviorLabel.LC_pair2 = LC;


        else
            behaveData_pair1 = behaveData(:,pair1StartInd:pair1EndInd);
            trialInfo_pair1 = taskPerformance_onebox(behaveData_pair1,IO_Fs,sessInfo(sessInd),'pair1',p);
            trialInfo_pair1.onsetInd = trialInfo_pair1.onsetInd + pair1StartInd -1;
            trialInfo_pair1.offsetInd = trialInfo_pair1.offsetInd + pair1StartInd -1;
            trialInfo_pair1.onsetTs = IO_ts(trialInfo_pair1.onsetInd);
            trialInfo_pair1.offsetTs = IO_ts(trialInfo_pair1.offsetInd);

            matureTrialAll_pair1 = sum(trialInfo_pair1.mature == 1);
            matureTrialA_pair1 = sum(trialInfo_pair1.mature == 1 & trialInfo_pair1.odor == 1);
            matureTrialB_pair1 = sum(trialInfo_pair1.mature == 1 & trialInfo_pair1.odor == 2);
            matureTrialC_pair1 = sum(trialInfo_pair1.mature == 1 & trialInfo_pair1.odor == 3);
            matureTrialD_pair1 = sum(trialInfo_pair1.mature == 1 & trialInfo_pair1.odor == 4);

            % initiate manual water trial label
            trialInfo_pair1.manual = nan(trialInfo_pair1.num,1);  
            trialInfo_pair1.manual(trialInfo_pair1.mature == 1) = 0;
            % trialInfo_pair2.manual = nan(trialInfo_pair2.num,1);
            % trialInfo_pair2.manual(trialInfo_pair2.mature == 1) = 0;

            % assign manual and compare trial detection from labview event log
            % compare the two event streams (onebox and labview log)
            % event log can be both longer or shorter than onebox pair
            % length, we sometimes stop labview or cut pair
            % read event from eventLOG.txt file
            log_folder_pair1 = dir(fullfile(sessInfo(sessInd).npx_path, '*pair1'));
            % log_folder_pair2 = dir(fullfile(sessInfo(sessInd).npx_path, '*pair2'));
            if ~isempty(log_folder_pair1)
                logFile_pair1 = dir(fullfile(log_folder_pair1.folder,log_folder_pair1.name,'*eventLOG.txt'));
                log_pair1 = load(fullfile(logFile_pair1.folder,logFile_pair1.name));
                [evnt_pair1,results_pair1] = gngGetTS_check_error_AB_ABCD_labView(log_pair1);
                matureInd = trialInfo_pair1.mature == 1;


                % find matching trials
                match = findContinuousOdorMatch(trialInfo_pair1, evnt_pair1);
                matchPct = match.bestLen / sum(matureInd) * 100;
                fprintf('pair1 %2.1f trials found in lab view eventLog\n',matchPct);
                
                oneBoxOutcome = trialInfo_pair1.outcome(matureInd);
                oneBoxOutcome2 = oneBoxOutcome(match.oneboxStart:match.oneboxEnd);
                labViewOutcome = evnt_pair1.action(match.eventStart:match.eventEnd);
                oneBoxOutcome2 = oneBoxOutcome2(:);
                labViewOutcome = labViewOutcome(:);

                oneBoxManual = zeros(sum(matureInd),1);
                oneBoxManual(match.oneboxStart:match.oneboxEnd) = evnt_pair1.MW(match.eventStart:match.eventEnd);
                trialInfo_pair1.manual(matureInd) = oneBoxManual;

                % for manual water trials, always assign outcome as 2
                % or 4, dependeing on odors
                % oneBoxOdor = trialInfo_pair2.odor(matureInd);
                % oneBoxOdor2 = oneBoxOdor(match.oneboxStart:match.oneboxEnd);
                % manualInd = oneBoxManual == 1;
                % odor 1 or 3 + manual water = outcome 2
                % oneBoxOutcome2(manualInd & ismember(oneBoxOdor2, [1 3])) = 2;
                % % odor 2 or 4 + manual water = outcome 3
                % oneBoxOutcome2(manualInd & ismember(oneBoxOdor2, [2 4])) = 3;

                % Assign back
                % oneBoxOutcome(match.oneboxStart:match.oneboxEnd) = oneBoxOutcome2;
                % trialInfo_pair2.outcome(matureInd) = oneBoxOutcome;
                % 
                % if any(oneBoxOutcome2(:) ~= labViewOutcome(:))
                %     error('Pair 2 Outcome mismatch between OneBox and LabView eventLog.');
                % else
                %     fprintf('pair2 %2.1f trials found in lab view eventLog\n',matchPct);
                %     fprintf('pair2 all trial outcome matches between onebox and labview\n');
                % end

            else
                warning('No labview_pair1 event folder')
            end

            behaviorLabel.trialInfo_pair1 = trialInfo_pair1;
            behaviorLabel.trialInfo_pair2 = [];

            h = figure(1);
            h.Position = [100,100,1200,900];
            subplot(4,1,1)
            f(1) = subplot(4,1,1);
            ts = IO_ts(pair1StartInd:pair1EndInd);
            plot(ts,behaveData_pair1(sessInfo(sessInd).odorA_ch,:))
            hold on
            plot(ts,behaveData_pair1(sessInfo(sessInd).odorB_ch,:))
            plot(ts,behaveData_pair1(sessInfo(sessInd).odorC_ch,:))
            plot(ts,behaveData_pair1(sessInfo(sessInd).odorD_ch,:))
            % plot(ts,data(sessInfo(sessInd).odorE_ch,:))
            % plot(ts,data(sessInfo(sessInd).odorF_ch,:))
            plot([ts(1),ts(end)],[p.odorThres,p.odorThres],'k')
            text(ts(1),p.odorThres+0.5,'Thres')
            legend({'A','B','C','D'})
            ylabel('Odor Amp (V)')
            TITLE1 = sprintf('%s-%s-pair1-onebox-IO-signal',behaviorLabel.animalID,behaviorLabel.recDate);
            title(TITLE1,'Interpreter','None');


            f(2) = subplot(4,1,2);
            plot(ts,behaveData_pair1(sessInfo(sessInd).led_ch,:))
            hold on
            plot([ts(1),ts(end)],[p.LEDThres,p.LEDThres],'k')
            text(ts(1),p.LEDThres+0.5,'Thres')
            ylabel('LED Amp (V)')
            ylim([0 4])

            f(3) = subplot(4,1,3);
            plot(ts,behaveData_pair1(sessInfo(sessInd).lick_ch,:))
            hold on
            plot([ts(1),ts(end)],[p.lickThres,p.lickThres],'k')
            text(ts(1),p.lickThres+0.5,'Thres')
            % label manual trials
            onsetInd = trialInfo_pair1.onsetInd;
            maualInd = onsetInd(trialInfo_pair1.manual == 1);
            if ~isempty(maualInd)
                plot(IO_ts(maualInd), 4.5*ones(size(ts(maualInd))),'r*')
            end
            ylabel('Lick Amp (V)')
            title('Lick signal not used in deciding correct or not, manual water *')
            ylim([0 5])

            f(4) = subplot(4,1,4);
            plot(ts,behaveData_pair1(sessInfo(sessInd).sucrose_ch,:))
            hold on
            plot(ts,behaveData_pair1(sessInfo(sessInd).quinine_ch,:))
            plot([ts(1),ts(end)],[p.rewardThres,p.rewardThres],'k')
            text(ts(1),p.rewardThres+0.5,'Thres')
            legend({'S','Q'})
            ylabel('Reward Amp (V)')
            xlabel('Time (s)')

            linkaxes(f,'x')
            figName = sprintf('%s%s%s-%s%s',savedir_plot,'\',behaviorLabel.animalID,behaviorLabel.recDate,'-pair1-onebox-IO-signal');
            print(figName,'-dpng','-r200');


            % plot the accuracy plot
            % plot the learning curve for pair 1
            h = figure(3);
            h.Position = [100,100,1200,900];
            subplot(2,1,1)
            [LC,matureRate] = learningCurveCalc(behaviorLabel.trialInfo_pair1,p.LC_window,p.LC_window_each);

            subplot(2,2,1)
            % AB vs CD
            plot(LC.ABCD(1,:),'b','LineWidth',2);hold on;axis tight
            plot(LC.ABCD(4,:),'r','LineWidth',2)
            % plot(matureRate.all,'k','LineWidth',2)
            ylim([30 100])
            label = ['Trials (mov window = ' num2str(p.LC_window) ' trials)'];
            xlabel(label)
            ylabel('Correct %')
            legend('Correct% A/B ','Correct% novel','Location','southeast')
            TITLE1 = sprintf('%s-%s-%s',behaviorLabel.animalID,behaviorLabel.recDate,cell2mat(behaviorLabel.odorPair1));
            TITLE2 = sprintf('Success trial num: %d',matureTrialAll_pair1);
            title({TITLE1;TITLE2},'interpreter','None')


            subplot(2,2,2)
            % A, B, C, D in all type trials
            plot(LC.ABCD(5,:),'Color',cTablegreen,'LineWidth',2);hold on;axis tight
            plot(LC.ABCD(6,:),'Color',cTableOrange,'LineWidth',2)
            plot(LC.ABCD(2,:),'r','LineWidth',2)
            plot(LC.ABCD(3,:),'b','LineWidth',2)
            ylim([0 100])
            label = ['Trials (mov window = ' num2str(p.LC_window) ' trials)'];
            xlabel(label)
            ylabel('Correct %')
            title('Each odor')
            legend('1','2','A','B','Location','southeast')

            subplot(2,2,3)
            plot(LC.C,'Color',cTablegreen,'LineWidth',2);
            hold on
            axis tight
            plot(LC.D,'Color',cTableOrange,'LineWidth',2);
            plot(LC.A,'r','LineWidth',2);
            plot(LC.B,'b','LineWidth',2);
            label = ['Trials (each odor); movWin = ' num2str(p.LC_window_each) ' trials'];
            xlabel(label)
            legend('1','2','A','B','Location','southeast')
            ylim([0 100])
            TITLE1 = sprintf('Trial num A: %d B: %d',matureTrialA_pair1,matureTrialB_pair1);
            TITLE2 = sprintf('Trial num 1: %d 2: %d',matureTrialC_pair1,matureTrialD_pair1);
            title({TITLE1;TITLE2},'interpreter','None')
            
            subplot(2,2,4)
            plot(LC.CD,'r','LineWidth',2);hold on;axis tight
            plot(LC.AB,'b','LineWidth',2);
            label = ['Trials (each odor); movWin = ' num2str(p.LC_window_each) ' trials'];
            xlabel(label)
            legend('12','AB','Location','southeast')
            ylim([0 100])
            TITLE1 = sprintf('Trial num AB %d',matureTrialA_pair1+matureTrialB_pair1);
            TITLE2 = sprintf('Trial num 1-2: %d',matureTrialC_pair1+matureTrialD_pair1);
            title({TITLE1;TITLE2},'interpreter','None')

            figName = sprintf('%s%s%s-%s%s',savedir_plot,'\',behaviorLabel.animalID,behaviorLabel.recDate,'-pair1-learningCurve');
            print(figName,'-dpng','-r200');

            behaviorLabel.LC_pair1 = LC;
            behaviorLabel.LC_pair2 = [];
        end
    end
    

    if p.saveFile
        save(fullfile(savedir,'behaviorLabel.mat'), 'behaviorLabel','-v7.3');
        % save(fullfile(sessInfo(sessInd).npx_path,'behaviorLabel.mat'), 'behaviorLabel','-v7.3');
    end
    fprintf('Finish behavioral analysis: %d\n', sessInd);
    clear behaviorLabel

end
close all
end
