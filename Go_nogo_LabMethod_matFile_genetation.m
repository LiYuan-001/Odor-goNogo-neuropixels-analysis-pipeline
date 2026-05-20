% this code is to save the .mats required when using Igarashi lab original
% code, for neuropixels analysis
% I tried keep original codes as much as possible
% Results here are used or code labeled with _LabMethod, not with other
% pipeline
% I did not organize this code well. This is just for temporary use
% Li YUAN, Tohoku, 2026-May

function Go_nogo_LabMethod_matFile_genetation(inFile,analyzeSes)
close all
% set(0, 'DefaultFigureVisible', 'off');
p.saveFile = 1;
p.savePlot = 1;
p.saveDir = 'N:\yuan\AnalysisResults'; % change it to your own folder

p.bin = 50/10^3; % unit sec, time bin for calculate rate
premsec = 1; % unit sec, time shown before onset of odor
postmsec = 4;
p.odor_len = 1;
opt = 0;
% p.postTime = 4; % unit sec, time shown after LED + reward period

% Read in input information
sessInfo = SessInfoImport_Npx(inFile);

for sessInd = analyzeSes(1:end)

    close all

    % used for saving file name
    mouse = sessInfo(sessInd).animalID;
    owner = 'npx';
    mouseline = 'WT';
    date = sessInfo(sessInd).recDate;
    inhStimTag = 0;
    blockID_opt = [];
    blockID_inhStim = [];
    Area = 'HPC';
    dirnRec = [];
    numOdor.ABsession = 2;
    numOdor.ABCDsession = 4;

    p.delay_len = sessInfo(sessInd).delay_len;
    p.response_len = sessInfo(sessInd).response_len;

    if p.saveFile
        savedir = sprintf('%s%s%s%s%s',p.saveDir,'\ProcessedData\',sessInfo(sessInd).animalID,'-date-',sessInfo(sessInd).recDate);
        if ~exist(savedir, 'dir')
            mkdir(savedir);
        end
        savedir2 = sprintf('%s%s%s%s%s%s',p.saveDir,'\ProcessedData\',sessInfo(sessInd).animalID,'-date-',sessInfo(sessInd).recDate,'\LabMethod');
        if ~exist(savedir2, 'dir')
            mkdir(savedir2);
        end
        savedir3 = sprintf('%s%s%s%s%s%s',p.saveDir,'\ProcessedData\',sessInfo(sessInd).animalID,'-date-',sessInfo(sessInd).recDate,'\LabMethod\psth');
        if exist(savedir3, 'dir')
            delete(fullfile(savedir3, '*'));
        else
            mkdir(savedir3);
        end
        savedir4 = sprintf('%s%s%s%s%s%s',p.saveDir,'\ProcessedData\',sessInfo(sessInd).animalID,'-date-',sessInfo(sessInd).recDate,'\LabMethod\spikeFile');
        if exist(savedir4, 'dir')
            delete(fullfile(savedir4, '*'));
        else
            mkdir(savedir4);
        end

    end

    
    % % load corrected TIMESTAMPS
    % syncFile = fullfile(savedir,'syncTime.mat');
    % syncTime = load(syncFile);
    % IO_ts = syncTime.Npx_timeStamps.IO.Ts_Sync;

    % get behavior matrix same as lab code use
    % from gngGetTS_ss_check_error_AB_ABCD
    % all trials
    % early trials
    % [1]odor type ,[2]odor onset, [3]NaN, [4]NaN
    % tmp_ss = [events(n_start + 3) events(n_start + 2) NaN NaN];
    % signal_ss=[signal_ss; tmp_ss];
    % mature trials
    % [1]odor type ,[2]odor onset, [3]Correct / Error, [4]manual_rwd
    % tmp_ss = [events(n_start + 3) events(n_start + 2) CoEr evnt.MW(i);];
    % signal_ss=[signal_ss; tmp_ss];
    % from gngGetTS2
    % Nlx_t_type_AB=Nlx_t_type; %  [time + odorID] in Nlx event file
    % TS_AB_Nlx=t;        %Nlx timestamp for mature odor onset
    % TS_id_AB=t_id;      %Nlx timestamp index for mature odor onset
    % TS_id_odor_AB=t_odor; %odor ID for mature odor onset (1-4 for A, B, C, D)

    % read in behavior time
    behaveFile = fullfile(savedir,'behaviorLabel.mat');
    load(behaveFile); % output: behaviorLabel
    % assuming always pair1 is 2 odor and pair2 is 4 odor
    TS_AB = [];
    TS_ABCD = [];
    Nlx_t_type_AB = [];
    Nlx_t_type_ABCD = [];
    for pairInd = 1:2 % assuming maximum 2 pairs each day
        pairName = sprintf('pair%d',pairInd);
        trialInfoName = sprintf('trialInfo_pair%d',pairInd);

        if ~isempty(sessInfo(sessInd).(pairName))

            odorName = sessInfo(sessInd).(pairName);
            odorAll = behaviorLabel.(trialInfoName).odor; % 1,2,3,4
            odorNum = length(unique(odorAll));
            onsetIndAll = behaviorLabel.(trialInfoName).onsetInd;
            offsetIndAll = behaviorLabel.(trialInfoName).offsetInd;
            onset_TsAll = behaviorLabel.(trialInfoName).onsetTs;
            offset_TsAll = behaviorLabel.(trialInfoName).offsetTs;
            outcomeLabelAll = behaviorLabel.(trialInfoName).outcome;
            manualLabelAll = (behaviorLabel.(trialInfoName).manual==1) & (behaviorLabel.(trialInfoName).mature==1);

            matureLabel = behaviorLabel.(trialInfoName).mature == 1; % 0, premature lick, 1, matureful trial
            outcomeLabel = outcomeLabelAll(matureLabel); % 1: Hit, 2: Miss, 3: FalseAlarm, 4: CorrectRejection
            odor = odorAll(matureLabel);
            onsetInd = onsetIndAll(matureLabel);
            offsetInd = offsetIndAll(matureLabel);
            onset_Ts = onset_TsAll(matureLabel);
            offset_Ts = offset_TsAll(matureLabel);
            outcome = outcomeLabelAll(matureLabel);
            manualLabel = manualLabelAll(matureLabel);

            if odorNum == 2
                blockID_AB = pairInd;
                Nlx_t_type_AB = [onset_TsAll(:),odorAll(:)]; %  [time + odorID] in Nlx event file
                % TS_AB_Nlx = onset_Ts(:);        %Nlx timestamp for mature odor onset
                % TS_id_AB = onsetInd(:);      %Nlx timestamp index for mature odor onset
                % TS_id_odor_AB= odor(:); %odor ID for mature odor onset (1-4 for A, B, C, D)
                trialInd = 1:length(odorAll);
                signal_ss = [trialInd(:),odorAll(:),onset_TsAll(:),outcomeLabelAll(:),manualLabelAll(:)];
                TS_AB = signal_ss;

                [Nlx_TS,Nlx_O,P_task] = cal_TS_TP_onebox(Nlx_t_type_AB,signal_ss);
                Nlx_TS.OptFlag = NaN;
                Nlx_TS_AB = Nlx_TS;
                Nlx_O_AB = Nlx_O;
                P_task_AB = P_task;

            elseif odorNum == 4
                blockID_ABCD = pairInd;
                Nlx_t_type_ABCD = [onset_TsAll(:),odorAll(:)]; %  [time + odorID] in Nlx event file
                % TS_ABCD_Nlx = onset_Ts(:);        %Nlx timestamp for mature odor onset
                % TS_id_ABCD = onsetInd(:);      %Nlx timestamp index for mature odor onset
                % TS_id_odor_ABCD= odor(:); %odor ID for mature odor onset (1-4 for A, B, C, D)
                trialInd = 1:length(odorAll);
                signal_ss = [trialInd(:),odorAll(:),onset_TsAll(:),outcomeLabelAll(:),manualLabelAll(:)];
                TS_ABCD = signal_ss;

                [Nlx_TS,Nlx_O,P_task] = cal_TS_TP_onebox(Nlx_t_type_ABCD,signal_ss);
                Nlx_TS.OptFlag = NaN;
                Nlx_TS_ABCD = Nlx_TS;
                Nlx_O_ABCD = Nlx_O;
                P_task_ABCD = P_task;

                LC = behaviorLabel.LC_pair2.ABCD;
                LC_A = behaviorLabel.LC_pair2.A;
                LC_B = behaviorLabel.LC_pair2.B;
                LC_C = behaviorLabel.LC_pair2.C;
                LC_D = behaviorLabel.LC_pair2.D;
                LC_AB = behaviorLabel.LC_pair2.AB;
                LC_CD = behaviorLabel.LC_pair2.CD;
                s_LC_AB = LC_AB;
                s_LC_CD = LC_CD;
                changePt = findchangepts(LC_CD);

                behavResult.numGo_1   = sum(odor == 1);
                behavResult.numNoGo_1 = sum(odor == 2);
                behavResult.numGo_2   = sum(odor == 3);
                behavResult.numNoGo_2 = sum(odor == 4);

                behavResult.numHit_1 = sum((odor == 1) & (outcome == 1));
                behavResult.numCR_1  = sum((odor == 2) & (outcome == 4));
                behavResult.numHit_2 = sum((odor == 3) & (outcome == 1));
                behavResult.numCR_2  = sum((odor == 4) & (outcome == 4));

                behavResult.CR1=((behavResult.numHit_1 + behavResult.numCR_1)/(behavResult.numGo_1 + behavResult.numNoGo_1)) * 100;% Correct rate Familiar
                behavResult.CR2=((behavResult.numHit_2 + behavResult.numCR_2)/(behavResult.numGo_2 + behavResult.numNoGo_2)) * 100; % Correct rate Novel

               
                folderID = '_2';
                odorpair1 = sessInfo(sessInd).pair1;
                odorpair2 = sessInfo(sessInd).pair2;
                inhStimTag = [];
                dirnRec =[];
                dirnAB = [];
                dirnABCD = [];
                behavFileName=fullfile(savedir,'LabMethod',strcat('behavFile_', mouseline,'_mouse',mouse,'_',date,folderID,'_4odors.mat'));

                save(behavFileName, 'owner', 'mouse', 'date', 'folderID', 'odorpair2', 'mouseline', 'inhStimTag','dirnRec','dirnAB','dirnABCD',...
                    'behavResult','LC','mouseline','LC_A','LC_B','LC_C','LC_D','LC_AB','LC_CD','s_LC_AB','s_LC_CD','changePt','behavFileName', 'Area')

            else
                error('Odor num error')
            end
        end
    end

    

    % from Make_TimeStamp_list.m
    % history time stamp
    % [1] connect AB and ABCD data
    tmpEventID=[TS_AB;TS_ABCD]; %LabVIEW TS
    tmpTS1=Nlx_t_type_AB(:,1);% %Nlx TS in sec
    tmpTS2=Nlx_t_type_ABCD(:,1);% sec
    tmp1=Nlx_t_type_AB(:,2);
    tmp2=Nlx_t_type_ABCD(:,2);
    tmpTS1=[tmpTS1 tmp1];
    tmpTS2=[tmpTS2 tmp2];
    tmpTS=[tmpTS1;tmpTS2];
    % [2] detect the first novel odor presentation
    first_Novel_trial=find(tmpEventID(:,2)==3 | tmpEventID(:,2)==4,1,'first');

    % [3] make new time stamp matrix
    history_trial_id=[-first_Novel_trial+1:0 1:size(tmpEventID,1)-first_Novel_trial]; %Centers trial numbers so that 0 = first novel
    history.time_stamp=[history_trial_id(:) tmpTS tmpEventID]; %One row per trial; first column is history_trial ID's where 0 is initial novel

    history.time_stamp_name(1,1)={'Trial number triggered by initial novel odor'};
    history.time_stamp_name(1,2)={'[0] is initial novel odor'};
    history.time_stamp_name(2,1)={'Time stamp for odor onset(s)'};
    history.time_stamp_name(3,1)={'Odor type obtained from Nlx data'};
    history.time_stamp_name(4,1)={'Trial number '};
    history.time_stamp_name(5,1)={'Odor type obtained from Labview'};
    history.time_stamp_name(6,1)={'Time stamp for odor onset(s) from Labview'};
    history.time_stamp_name(7,1)={'ID for performance'};
    history.time_stamp_name(7,2)={'[1]Correct, [0]Incorrect, [NaN]immmature trial'};
    history.time_stamp_name(8,1)={'ID for manual rwd'};
    history.time_stamp_name(8,2)={'[0]No manual rwd, [1]Manual rwd, [NaN]immmature trial'};


    EventInfo = [];
    Nlx_TS=Nlx_TS_AB;
    Nlx_O=Nlx_O_AB;
    saveName=fullfile(savedir,'LabMethod',strcat('Nlx_TS_log_0',num2str(blockID_AB),'.mat')); % saving same file without _AB
    save(saveName, 'Nlx_TS','Nlx_O','EventInfo');
    Nlx_TS=Nlx_TS_ABCD;
    Nlx_O=Nlx_O_ABCD;
    saveName=fullfile(savedir,'LabMethod',strcat('Nlx_TS_log_0',num2str(blockID_ABCD),'.mat'));
    save(saveName, 'Nlx_TS','Nlx_O','EventInfo');

    f=dir(fullfile(savedir,'*_SessionInfo.mat'));
    if ~isempty(f)
        delete(fullfile(savedir,f.name)) %delete previous SessionInfo.mat files
    end

    mouse = sessInfo(sessInd).animalID;
    date = sessInfo(sessInd).recDate;
    odorpair1 = cell2mat(sessInfo(sessInd).pair1);
    odorpair2 = cell2mat(sessInfo(sessInd).pair2);
    behavFileName = [];
    sessionInfoFileName = fullfile(savedir,'LabMethod',strcat('sessionInfo_', mouseline,'_mouse',mouse,'_',date,'_',odorpair2,'.mat'));

    save(sessionInfoFileName, 'owner', 'mouse', 'date', 'odorpair1', 'odorpair2', 'mouseline', 'inhStimTag', 'blockID_opt',...   % these parameters are also saved in eventID.mat
        'blockID_AB','blockID_inhStim','blockID_ABCD','numOdor', 'Area',...                                                      % these parameters are also saved in eventID.mat
        'P_task_AB','P_task_ABCD','Nlx_t_type_AB','Nlx_t_type_ABCD','TS_AB','TS_ABCD',...                                        % added specifically in sessionInfo
        'dirnRec','history', 'behavFileName','sessionInfoFileName');                               % added specifically in sessionInfo

   % from eventFileChecker
   eventIDFileName = fullfile(savedir,'LabMethod','eventID.mat');
   save(eventIDFileName,'owner', 'mouse', 'date', 'odorpair1', 'odorpair2', 'mouseline', 'inhStimTag', 'blockID_opt','blockID_AB','blockID_inhStim','blockID_ABCD','numOdor', 'Area')

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

        % read in cluster spike time
        spk_file = sprintf('%s%s%s',savedir,'\npx_Cluster_',probeName);
        load(spk_file); % output:npx_Cluster

        clusterNum = npx_Cluster.clusterNum;
        xcoord = nan(clusterNum,1);
        ycoord = nan(clusterNum,1);
        cluster_ksID = nan(clusterNum,1);

        % from gngRaster2Odor2
        pairInd = 1;
        pairName = sprintf('pair%d',pairInd);
        if ~isempty(sessInfo(sessInd).(pairName))

            odorType_AB{1}=sort([Nlx_TS_AB.Hit1;Nlx_TS_AB.Miss1]);
            odorType_AB{2}=sort([Nlx_TS_AB.CR2;Nlx_TS_AB.FA2]);
            for i = 1:clusterNum

                tmp_spike = npx_Cluster.clusterInfo{i}.spkTs;
                xcoord(i) = npx_Cluster.clusterInfo{i}.channel.xcoords;
                ycoord(i) = npx_Cluster.clusterInfo{i}.channel.ycoords;
                cluster_ksID(i) = npx_Cluster.clusterInfo{i}.ID;

                for ii = 1:2 %repeat across trial types
                    tmp_ev = odorType_AB{ii};
                    data_total =[];
                    FR=[];
                    tmpname = strcat('odor',sessInfo(sessInd).pair1(ii));

                    for k = 1:length(tmp_ev)
                        piyo = tmp_spike-(tmp_ev(k));
                        data_total(k).times = piyo((-premsec < piyo) & (piyo < postmsec))'*1000;
                        FR(k).ITI = length(piyo((-5 < piyo) & (piyo < 0)))/5; % ITI = 5sec
                        FR(k).Odor = length(piyo((0 < piyo) & (piyo < p.odor_len)))/1;
                        FR(k).delay = length(piyo((p.odor_len < piyo) & (piyo < (p.odor_len+p.delay_len))))/2;
                    end

                    trialNum=length(data_total);
                    edges=[-premsec:p.bin:postmsec]; %Defines time windows for calculating rate
                    numSpikes=zeros(trialNum,round(1+(premsec+postmsec)/p.bin));
                    for k=1:trialNum
                        numSpikes(k,:)=histc(data_total(k).times,edges);
                    end
                    binRate = numSpikes/(p.bin);

                    if trialNum>10
                        numInitTrials=10; %
                        numLastTrials=10;
                    else
                        numInitTrials=trialNum; %
                        numLastTrials=trialNum;
                    end

                    binRateInitial = numSpikes(1:numInitTrials,:)/p.bin;
                    binRateEnd = numSpikes(end-numLastTrials+1:end,:)/p.bin;
                    binRateAll = numSpikes(1:end,:)/p.bin;

                    rate2(1,:)=mean(binRateInitial,1);
                    rate2(2,:)=mean(binRateEnd,1);
                    rate2(3,:)=mean(binRateAll,1);

                    rateSE(1,:)=std(binRateInitial)./sqrt(numInitTrials);
                    rateSE(2,:)=std(binRateEnd)./sqrt(numLastTrials);
                    rateSE(3,:)=std(binRateAll)./sqrt(size(numSpikes,1));

                    rate2S2{ii} =smooth_gaussian2(rate2,100,p.bin*10^3,1);
                    rateSE2{ii} =smooth_gaussian2(rateSE,100,p.bin*10^3,1);

                    PSTH_avg(ii) = {smooth_gaussian2(rate2,100,p.bin*10^3,1)};
                    PSTH_se(ii) = {smooth_gaussian2(rateSE,100,p.bin*10^3,1)};
                    raster(ii) = {data_total};

                end
                task_time.Nlx_O=Nlx_O_ABCD;
                task_time.Nlx_TS=Nlx_TS_ABCD;
                posAxisPSTH = [];

                data_name.raster={'[1]Go_Familiar, [2]Nogo_Familiar, [3]Go_Novel, [4]Nogo_Novel'};
                data_name.PSTH10trials={'[1]First_10trials, [2]Last_10trials,[3]All'};
                data_name.PSTH5trials={'[1]First_5trials, [2]Last_5trials,[3]All'};

                outf_PSTH_AB = sprintf('%s%s%s%s%s%s%d%s',savedir,'\LabMethod\psth\',sessInfo(sessInd).animalID,'_',sessInfo(sessInd).recDate,'_cluster_',i,'_PSTH_AB.mat');
                %Save TT.PSTH file: note this is for ABCD session in the second block. TODO:rename with a clearer name
                save(outf_PSTH_AB,'raster','posAxisPSTH','PSTH_avg','PSTH_se','data_name','task_time');

            end
        end

        % from gngRaster4Odor2
        pairInd = 2;
        pairName = sprintf('pair%d',pairInd);
        if ~isempty(sessInfo(sessInd).(pairName))
            odorType_ABCD{1}=sort([Nlx_TS_ABCD.Hit1;Nlx_TS_ABCD.Miss1]);
            odorType_ABCD{2}=sort([Nlx_TS_ABCD.CR2;Nlx_TS_ABCD.FA2]);
            odorType_ABCD{3}=sort([Nlx_TS_ABCD.Hit3;Nlx_TS_ABCD.Miss3]);
            odorType_ABCD{4}=sort([Nlx_TS_ABCD.CR4;Nlx_TS_ABCD.FA4]);

            for i = 1:clusterNum

                tmp_spike = npx_Cluster.clusterInfo{i}.spkTs;
                xcoord(i) = npx_Cluster.clusterInfo{i}.channel.xcoords;
                ycoord(i) = npx_Cluster.clusterInfo{i}.channel.ycoords;
                cluster_ksID(i) = npx_Cluster.clusterInfo{i}.ID;

                for ii = 1:4 %repeat across trial types
                    tmp_ev = odorType_ABCD{ii};
                    data_total =[];
                    FR=[];
                    tmpname = strcat('odor',sessInfo(sessInd).pair2(ii));

                    for k = 1:length(tmp_ev)
                        piyo = tmp_spike-(tmp_ev(k));
                        data_total(k).times = piyo((-premsec < piyo) & (piyo < postmsec))'*1000;
                        FR(k).ITI = length(piyo((-5 < piyo) & (piyo < 0)))/5; % ITI = 5sec
                        FR(k).Odor = length(piyo((0 < piyo) & (piyo < p.odor_len)))/1;
                        FR(k).delay = length(piyo((p.odor_len < piyo) & (piyo < (p.odor_len+p.delay_len))))/2;
                    end

                    trialNum=length(data_total);
                    edges=[-premsec:p.bin:postmsec]; %Defines time windows for calculating rate
                    numSpikes=zeros(trialNum,round(1+(premsec+postmsec)/p.bin));
                    for k=1:trialNum
                        numSpikes(k,:)=histc(data_total(k).times,edges);
                    end
                    binRate = numSpikes/(p.bin);

                    if trialNum>10
                        numInitTrials=10; %
                        numLastTrials=10;
                    else
                        numInitTrials=trialNum; %
                        numLastTrials=trialNum;
                    end

                    binRateInitial = numSpikes(1:numInitTrials,:)/p.bin;
                    binRateEnd = numSpikes(end-numLastTrials+1:end,:)/p.bin;
                    binRateAll = numSpikes(1:end,:)/p.bin;

                    rate2(1,:)=mean(binRateInitial,1);
                    rate2(2,:)=mean(binRateEnd,1);
                    rate2(3,:)=mean(binRateAll,1);

                    rateSE(1,:)=std(binRateInitial)./sqrt(numInitTrials);
                    rateSE(2,:)=std(binRateEnd)./sqrt(numLastTrials);
                    rateSE(3,:)=std(binRateAll)./sqrt(size(numSpikes,1));

                    PSTH_avg(ii) = {smooth_gaussian2(rate2,100,p.bin*10^3,1)};
                    PSTH_se(ii) = {smooth_gaussian2(rateSE,100,p.bin*10^3,1)};
                    raster(ii) = {data_total};

                end
                task_time.Nlx_O=Nlx_O_ABCD;
                task_time.Nlx_TS=Nlx_TS_ABCD;
                posAxisPSTH = [];

                data_name.raster={'[1]Go_Familiar, [2]Nogo_Familiar, [3]Go_Novel, [4]Nogo_Novel'};
                data_name.PSTH10trials={'[1]First_10trials, [2]Last_10trials,[3]All'};
                data_name.PSTH5trials={'[1]First_5trials, [2]Last_5trials,[3]All'};

                outf_PSTH_ABCD = sprintf('%s%s%s%s%s%s%d%s',savedir,'\LabMethod\psth\',sessInfo(sessInd).animalID,'_',sessInfo(sessInd).recDate,'_cluster_',i,'_PSTH.mat');
                %Save TT.PSTH file: note this is for ABCD session in the second block. TODO:rename with a clearer name
                save(outf_PSTH_ABCD,'raster','posAxisPSTH','PSTH_avg','PSTH_se','data_name','task_time');


                % from gng_cal_history_V5 & gngSessionSum3                
                % i am tired so not gonna organize code anymore -- Li YUAN
                numTrial2=20;
                numTria1=numTrial2+20;
                perfm_f=mean(LC(1,end-numTria1:end-numTrial2)); %Correct rate AB (see TaskPerformance_4odors)
                perfm_n=mean(LC(4,end-numTria1:end-numTrial2)); %Correct rate CD

                if perfm_f>=75 && perfm_n>=75
                    learning_idx=0;
                    % disp('Good performance')
                elseif perfm_f>=75 && perfm_n<75% familiar good
                    learning_idx=1;
                    % disp('Bad performance for novel pairs')
                elseif perfm_f<75 && perfm_n>=75% novelty good
                    learning_idx=2;
                    % disp('Bad performance for familiar pairs')
                else
                    learning_idx=3;
                    % disp('Bad performance')
                end

                numTrial = size(LC, 2);
                if numTrial >160
                    trialRange = 120:160; % 200815 Correct% will be calculated from trials 120-160
                else
                    trialRange = numTrial-39:numTrial;
                end

                perfm_f=mean(LC(1,trialRange));
                perfm_n=mean(LC(4,trialRange));

                learning_idx2(1) = perfm_f; % correct% for AB
                learning_idx2(2) = perfm_n; % correct% for CD

                TS.odorA=sort([Nlx_TS_ABCD.Hit1;Nlx_TS_ABCD.Miss1]);
                TS.odorB=sort([Nlx_TS_ABCD.CR2;Nlx_TS_ABCD.FA2]);
                TS.odorC=sort([Nlx_TS_ABCD.Hit3;Nlx_TS_ABCD.Miss3]);
                TS.odorD=sort([Nlx_TS_ABCD.CR4;Nlx_TS_ABCD.FA4]);

                outf_PSTH_AB = sprintf('%s%s%s%s%s%s%d%s',savedir,'\LabMethod\psth\',sessInfo(sessInd).animalID,'_',sessInfo(sessInd).recDate,'_cluster_',i,'_PSTH_AB.mat');              
                load(outf_PSTH_AB)
                ABdata=load(outf_PSTH_ABCD);
                if opt==1
                    load(f_opt(ss).name)
                end

                for yy=1:size(raster,2)
                    tmprast=raster{yy};
                    tmp_rast=[];
                    for ys=1:length(tmprast)
                        tmpd=tmprast(ys).times;
                        tmp=[ones(length(tmpd),1)*ys tmpd'];
                        tmp_rast=[tmp_rast;tmp];
                    end
                    raster2{yy}=tmp_rast;
                end

                for yy=1:size(ABdata.raster,2)
                    tmprast=ABdata.raster{yy};
                    tmp_rast=[];
                    for ys=1:length(tmprast)
                        tmpd=tmprast(ys).times;
                        tmp=[ones(length(tmpd),1)*ys tmpd'];
                        tmp_rast=[tmp_rast;tmp];
                    end
                    ABdata.raster2{yy}=tmp_rast;
                end
               
                spike_time = tmp_spike;
                % numOdor.ABCDsession = 4;
                % tired. not gonna organize it anymore. 
                gng_cal_history_V5

                spikeW = [];
                halfAmpdur = [];
                meanFR = [];
                acg = [];
                isi = [];
                waveAve = [];
                PSTH_z = [];
                spikeFileName=sprintf('%s%s%s%s%s%s%s%s%s%s%s%s%d%s',savedir,'\LabMethod\spikeFile\spikeFile_',Area,'_',mouseline,'_',mouse,'_',date,'_',odorpair2,'_cluster_',i,'.mat');
                save(spikeFileName, 'owner', 'mouse', 'date', 'odorpair1', 'odorpair2', 'mouseline', 'inhStimTag', 'blockID_opt',...            % saved also in sessionInfo
                    'blockID_AB','blockID_inhStim','blockID_ABCD','numOdor', 'Area',...                                                         % saved also in sessionInfo
                    'P_task_AB','P_task_ABCD','Nlx_t_type_AB','Nlx_t_type_ABCD','TS_AB','TS_ABCD',...                                           % saved also in sessionInfo
                    'dirnRec','history','dirnRec','dirnAB','dirnABCD', 'behavFileName','sessionInfoFileName',...                                % saved also in sessionInfo
                    'behavResult','LC','mouseline','LC_A','LC_B','LC_C','LC_D','LC_AB','LC_CD','s_LC_AB','s_LC_CD','changePt',...               % saved also in behavFile
                    'spikeFileName','data_name','posAxisPSTH','PSTH_avg','PSTH_se','raster','spikeW','spike_time','halfAmpdur','meanFR',...     % saved specifically in spikeFile
                    'acg','isi','task_time','raster2','PSTH_z','TS','learning_idx','learning_idx2','waveAve','ABdata','opt')                    % saved specifically in spikeFile


            end


        end

    end

    fprintf('Finish lab original method matrix generation: %d\n', sessInd);

end
close all
end

function [smoothed2]=smooth_gaussian2(data,sigma,bin_width,dimension)
%sigma bin_width in sec (.004 for 4ms, bin_width=.001 for 1000Hz)

% dimension==1: smooth across columns for data in each rows
% dimension==2: smooth across raws for data in each columns

ratio=10*sigma/bin_width;


if dimension==1
    
    for i=1:size(data,1)
        data2=[repmat(data(i,1),1,ratio), data(i,:), repmat(data(i,size(data,2)),1,ratio)]; % prolong edges
        edges=[-3*sigma:bin_width:3*sigma];
        kernel=normpdf(edges,0,sigma);
        kernel=kernel*bin_width; %multiply by bin width
        center=ceil(length(edges)/2);
        smoothed=conv(data2,kernel);
        smoothed=smoothed(center:numel(data2)+center-1);
        smoothed2(i,:)=smoothed(ratio+1:numel(data2)-ratio); % cut edges to proper length
        
    end
    
    
    
elseif  dimension==2
    
    for i=1:size(data,2)
        data2=[repmat(data(1,i),ratio,1); data(:,i); repmat(data(size(data,1),i),ratio,1)];
        edges=[-3*sigma:bin_width:3*sigma];
        kernel=normpdf(edges,0,sigma);
        kernel=kernel*bin_width; %multiply by bin width
        center=ceil(length(edges)/2);
        smoothed=conv(data2,kernel);
        smoothed=smoothed(center:numel(data2)+center-1);
        smoothed2(:,i)=smoothed(ratio+1:numel(data2)-ratio);
        
    end
    
    
end

end


