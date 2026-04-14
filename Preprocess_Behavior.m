% this code is to read in the IO recorded labview data and analyze behavior
% Li YUAN, 2026-Feb
function Preprocess_Behavior(inFile,analyzeSes)
close all
p.saveFile = 1;
p.savePlot = 1;
p.saveDir = 'N:\yuan\AnalysisResults'; % change it to your own folder

p.dataType = 'A';
p.LC_window = 20; % unit trials
p.LC_window_each = 5;

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

        IO_ts = (0:size(behaveData,2)-1) / IO_Fs;

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
            trialInfo_pair1 = taskPerformance_onebox(behaveData_pair1,IO_Fs,sessInfo(sessInd),'pair1');
            trialInfo_pair2 = taskPerformance_onebox(behaveData_pair2,IO_Fs,sessInfo(sessInd),'pair2');
            trialInfo_pair1.onsetInd = trialInfo_pair1.onsetInd + pair1StartInd -1;
            trialInfo_pair1.offsetInd = trialInfo_pair1.offsetInd + pair1StartInd -1;
            trialInfo_pair2.onsetInd = trialInfo_pair2.onsetInd + pair2StartInd -1;
            trialInfo_pair2.offsetInd = trialInfo_pair2.offsetInd + pair2StartInd -1;
            behaviorLabel.trialInfo_pair1 = trialInfo_pair1;
            behaviorLabel.trialInfo_pair2 = trialInfo_pair2;

            successTrialAll_pair1 = sum(behaviorLabel.trialInfo_pair1.success == 1);
            successTrialA_pair1 = sum(behaviorLabel.trialInfo_pair1.success == 1 & behaviorLabel.trialInfo_pair1.Odor == 1);
            successTrialB_pair1 = sum(behaviorLabel.trialInfo_pair1.success == 1 & behaviorLabel.trialInfo_pair1.Odor == 2);
            successTrialC_pair1 = sum(behaviorLabel.trialInfo_pair1.success == 1 & behaviorLabel.trialInfo_pair1.Odor == 3);
            successTrialD_pair1 = sum(behaviorLabel.trialInfo_pair1.success == 1 & behaviorLabel.trialInfo_pair1.Odor == 4);

            successTrialAll_pair2 = sum(behaviorLabel.trialInfo_pair2.success == 1);
            successTrialA_pair2 = sum(behaviorLabel.trialInfo_pair2.success == 1 & behaviorLabel.trialInfo_pair2.Odor == 1);
            successTrialB_pair2 = sum(behaviorLabel.trialInfo_pair2.success == 1 & behaviorLabel.trialInfo_pair2.Odor == 2);
            successTrialC_pair2 = sum(behaviorLabel.trialInfo_pair2.success == 1 & behaviorLabel.trialInfo_pair2.Odor == 3);
            successTrialD_pair2 = sum(behaviorLabel.trialInfo_pair2.success == 1 & behaviorLabel.trialInfo_pair2.Odor == 4);

            h = figure(1);
            h.Position = [100,100,1200,900];
            subplot(4,1,1)
            TITLE1 = sprintf('%s-%s-pair1-onebox-IO-signal',behaviorLabel.animalID,behaviorLabel.recDate);
            title(TITLE1,'Interpreter','None');
            figName = sprintf('%s%s%s-%s%s',savedir_plot,'\',behaviorLabel.animalID,behaviorLabel.recDate,'-pair1-onebox-IO-signal');
            print(figName,'-dpng','-r200');

            h = figure(2);
            h.Position = [100,100,1200,900];
            subplot(4,1,1)
            TITLE1 = sprintf('%s-%s-pair2-onebox-IO-signal',behaviorLabel.animalID,behaviorLabel.recDate);
            title(TITLE1,'Interpreter','None');
            figName = sprintf('%s%s%s-%s%s',savedir_plot,'\',behaviorLabel.animalID,behaviorLabel.recDate,'-pair2-onebox-IO-signal');
            print(figName,'-dpng','-r200');

            % plot the accuracy plot
            % plot the learning curve for pair 1
            h = figure(3);
            h.Position = [100,100,1200,900];
            subplot(2,1,1)
            [LC,successRate] = learningCurveCalc(behaviorLabel.trialInfo_pair1,p.LC_window,p.LC_window_each);

            subplot(2,2,1)
            % AB vs CD
            plot(LC.ABCD(1,:),'b','LineWidth',2);hold on;axis tight
            plot(LC.ABCD(4,:),'r','LineWidth',2)
            plot(successRate,'k','LineWidth',2)
            ylim([30 100])
            label = ['Trials (mov window = ' num2str(p.LC_window) ' trials)'];
            xlabel(label)
            ylabel('Correct %')
            legend('Correct% A/B ','Correct% novel','Success Rate','Location','southeast')
            TITLE1 = sprintf('%s-%s-%s',behaviorLabel.animalID,behaviorLabel.recDate,cell2mat(behaviorLabel.odorPair1));
            TITLE2 = sprintf('Success trial num: %d',successTrialAll_pair1);
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
            TITLE1 = sprintf('Trial num A: %d B: %d',successTrialA_pair1,successTrialB_pair1);
            TITLE2 = sprintf('Trial num 1: %d 2: %d',successTrialC_pair1,successTrialD_pair1);
            title({TITLE1;TITLE2},'interpreter','None')
            
            subplot(2,2,4)
            plot(LC.CD,'r','LineWidth',2);hold on;axis tight
            plot(LC.AB,'b','LineWidth',2);
            label = ['Trials (each odor); movWin = ' num2str(p.LC_window_each) ' trials'];
            xlabel(label)
            legend('12','AB','Location','southeast')
            ylim([0 100])
            TITLE1 = sprintf('Trial num AB %d',successTrialA_pair1+successTrialB_pair1);
            TITLE2 = sprintf('Trial num 1-2: %d',successTrialC_pair1+successTrialD_pair1);
            title({TITLE1;TITLE2},'interpreter','None')
            figName = sprintf('%s%s%s-%s%s',savedir_plot,'\',behaviorLabel.animalID,behaviorLabel.recDate,'-pair1-learningCurve');
            print(figName,'-dpng','-r200');

            behaviorLabel.LC_pair1 = LC;


            % plot the learning curve for pair 2
            h = figure(4);
            h.Position = [100,100,1200,900];
            subplot(2,1,1)
            [LC,successRate] = learningCurveCalc(behaviorLabel.trialInfo_pair2,p.LC_window,p.LC_window_each);

            subplot(2,2,1)
            % AB vs CD
            plot(LC.ABCD(1,:),'b','LineWidth',2);hold on;axis tight
            plot(LC.ABCD(4,:),'r','LineWidth',2)
            plot(successRate,'k','LineWidth',2)
            ylim([30 100])
            label = ['Trials (mov window = ' num2str(p.LC_window) ' trials)'];
            xlabel(label)
            ylabel('Correct %')
            legend('Correct% A/B ','Correct% novel','Success Rate','Location','southeast')
            TITLE1 = sprintf('%s-%s-%s',behaviorLabel.animalID,behaviorLabel.recDate,cell2mat(behaviorLabel.odorPair2));
            TITLE2 = sprintf('Success trial num: %d',successTrialAll_pair2);
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
            TITLE1 = sprintf('Trial num A: %d B: %d',successTrialA_pair2,successTrialB_pair2);
            TITLE2 = sprintf('Trial num 1: %d 2: %d',successTrialC_pair2,successTrialD_pair2);
            title({TITLE1;TITLE2},'interpreter','None')
            
            subplot(2,2,4)
            plot(LC.CD,'r','LineWidth',2);hold on;axis tight
            plot(LC.AB,'b','LineWidth',2);
            label = ['Trials (each odor); movWin = ' num2str(p.LC_window_each) ' trials'];
            xlabel(label)
            legend('12','AB','Location','southeast')
            ylim([0 100])
            TITLE1 = sprintf('Trial num AB %d',successTrialA_pair2+successTrialB_pair2);
            TITLE2 = sprintf('Trial num 1-2: %d',successTrialC_pair2+successTrialD_pair2);
            title({TITLE1;TITLE2},'interpreter','None')
            figName = sprintf('%s%s%s-%s%s',savedir_plot,'\',behaviorLabel.animalID,behaviorLabel.recDate,'-pair2-learningCurve');
            print(figName,'-dpng','-r200');

            behaviorLabel.LC_pair2 = LC;


        else
            behaveData_pair1 = behaveData(:,pair1StartInd:pair1EndInd);
            trialInfo_pair1 = taskPerformance_onebox(behaveData_pair1,IO_Fs,sessInfo(sessInd),'pair1');
            trialInfo_pair1.onsetInd = trialInfo_pair1.onsetInd + pair1StartInd -1;
            trialInfo_pair1.offsetInd = trialInfo_pair1.offsetInd + pair1StartInd -1;
            behaviorLabel.trialInfo_pair1 = trialInfo_pair1;
            behaviorLabel.trialInfo_pair2 = [];

            successTrialAll_pair1 = sum(behaviorLabel.trialInfo_pair1.success == 1);
            successTrialA_pair1 = sum(behaviorLabel.trialInfo_pair1.success == 1 & behaviorLabel.trialInfo_pair1.odor == 1);
            successTrialB_pair1 = sum(behaviorLabel.trialInfo_pair1.success == 1 & behaviorLabel.trialInfo_pair1.odor == 2);
            successTrialC_pair1 = sum(behaviorLabel.trialInfo_pair1.success == 1 & behaviorLabel.trialInfo_pair1.odor == 3);
            successTrialD_pair1 = sum(behaviorLabel.trialInfo_pair1.success == 1 & behaviorLabel.trialInfo_pair1.odor == 4);

            h = figure(1);
            h.Position = [100,100,1200,900];
            subplot(4,1,1)
            TITLE1 = sprintf('%s-%s-pair1-onebox-IO-signal',behaviorLabel.animalID,behaviorLabel.recDate);
            title(TITLE1,'Interpreter','None');
            figName = sprintf('%s%s%s-%s%s',savedir,'\',behaviorLabel.animalID,behaviorLabel.recDate,'-pair1-onebox-IO-signal');
            print(figName,'-dpng','-r200');

            % plot the accuracy plot
            % plot the learning curve for pair 1
            h = figure(3);
            h.Position = [100,100,1200,900];
            subplot(2,1,1)
            [LC,successRate] = learningCurveCalc(behaviorLabel.trialInfo_pair1,p.LC_window,p.LC_window_each);

            subplot(2,2,1)
            % AB vs CD
            plot(LC.ABCD(1,:),'b','LineWidth',2);hold on;axis tight
            plot(LC.ABCD(4,:),'r','LineWidth',2)
            % plot(successRate.all,'k','LineWidth',2)
            ylim([30 100])
            label = ['Trials (mov window = ' num2str(p.LC_window) ' trials)'];
            xlabel(label)
            ylabel('Correct %')
            legend('Correct% A/B ','Correct% novel','Success Rate','Location','southeast')
            TITLE1 = sprintf('%s-%s-%s',behaviorLabel.animalID,behaviorLabel.recDate,cell2mat(behaviorLabel.odorPair1));
            TITLE2 = sprintf('Success trial num: %d',successTrialAll_pair1);
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
            TITLE1 = sprintf('Trial num A: %d B: %d',successTrialA_pair1,successTrialB_pair1);
            TITLE2 = sprintf('Trial num 1: %d 2: %d',successTrialC_pair1,successTrialD_pair1);
            title({TITLE1;TITLE2},'interpreter','None')
            
            subplot(2,2,4)
            plot(LC.CD,'r','LineWidth',2);hold on;axis tight
            plot(LC.AB,'b','LineWidth',2);
            label = ['Trials (each odor); movWin = ' num2str(p.LC_window_each) ' trials'];
            xlabel(label)
            legend('12','AB','Location','southeast')
            ylim([0 100])
            TITLE1 = sprintf('Trial num AB %d',successTrialA_pair1+successTrialB_pair1);
            TITLE2 = sprintf('Trial num 1-2: %d',successTrialC_pair1+successTrialD_pair1);
            title({TITLE1;TITLE2},'interpreter','None')

            figName = sprintf('%s%s%s-%s%s',savedir_plot,'\',behaviorLabel.animalID,behaviorLabel.recDate,'-pair1-learningCurve');
            print(figName,'-dpng','-r200');

            behaviorLabel.LC_pair1 = LC;
            behaviorLabel.LC_pair2 = [];
        end
    end


    %% compare the two event streams (onebox and labview log) if both exist
    % read event from eventLOG.txt file if exist
    log_folder_pair1 = dir(fullfile(sessInfo(sessInd).npx_path, '*pair1'));
    log_folder_pair2 = dir(fullfile(sessInfo(sessInd).npx_path, '*pair2'));
    if ~isempty(log_folder_pair1)
        logFile_pair1 = dir(fullfile(fullfile(log_folder_pair1.folder,log_folder_pair1.name,'*eventLOG.txt')));
        log_pair1 = load(fullfile(logFile_pair1.folder,logFile_pair1.name));
        [evnt_pair1,results_pair1] = gngGetTS_check_error_AB_ABCD_labView(log_pair1);
        trialInfo_evnt = trialInfo_pair1.outcome';
        trialInfo_evnt = trialInfo_evnt(~isnan(trialInfo_evnt));
        trialNumMax = min(trialInfo_evnt,length(evnt_pair1.action));
        if any(trialInfo_pair1.success(1:trialNumMax)' ~= results_pair1.suc_delay(1:trialNumMax)) ||...
                any(trialInfo_evnt(1:trialNumMax) ~= evnt_pair1.action(1:trialNumMax))
            fprintf('pair1: onebox detection doesnt match log txt\n');
        else
            fprintf('pair1: onebox detection matches log txt\n');
        end
    end
    if ~isempty(log_folder_pair2)
        logFile_pair2 = dir(fullfile(log_folder_pair2.folder,log_folder_pair2.name,'*.eventLOG.txt'));
        log_pair2 = load(fullfile(logFile_pair2.folder,logFile_pair2.name));
        [evnt_pair2,results_pair2] = gngGetTS_check_error_AB_ABCD_labView(log_pair2);
        trialInfo_evnt = trialInfo_pair2.outcome';
        trialInfo_evnt = trialInfo_evnt(~isnan(trialInfo_evnt));
        trialNumMax = min(trialInfo_evnt,length(evnt_pair2.action));
        if any(trialInfo_pair2.success(1:trialNumMax)' ~= results_pair2.suc_delay(1:trialNumMax)) || ...
                any(trialInfo_evnt(1:trialNumMax) ~= evnt_pair2.action(1:trialNumMax))
            fprintf('pair2: onebox detection doesnt match log txt\n');
        else
            fprintf('pair2: onebox detection matches log txt\n');
        end
    end

    if p.saveFile
        save(fullfile(savedir,'behaviorLabel.mat'), 'behaviorLabel','-v7.3');
    end
    fprintf('Finish behavioral analysis: %d\n', sessInd);
    clear behaviorLabel

end
close all
end
