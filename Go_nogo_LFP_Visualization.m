% this code is to read in spike timestamp and align it to the go-nogo event
% time
% Li YUAN, Tohoku, 2026-Apr
function Go_nogo_LFP_Visualization(inFile,analyzeSes,LFP_ch)
close all
p.saveFile = 1;
p.savePlot = 1;
p.saveDir = 'N:\yuan\AnalysisResults'; % change it to your own folder

p.bin = 100/10^3; % unit sec, time bin for calculate rate
p.smoothTime = 100/10^3; % unit sec, time bin around central bin
p.preTime = 0; % unit sec, time shown before offset of odor
p.postTime = 0; % unit sec, time shown after LED + reward period
p.eventNum = 20;
p.plotStep = 600; % unit uV
LFP_ch = LFP_ch + 1;

p.space = 'lin'; % 'log' or 'lin'  spacing of f's
p.frange = [3 250];
p.nCycles = 10; % cycles for morelet wavelet
p.dt = 0.01; 
p.smoothTime = 0.05;
p.method = 'wavelet';
p.nfreqs = 100;

colorMap = [];

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

    goNogo_LFP_Matrix.animalID = sessInfo(sessInd).animalID;
    goNogo_LFP_Matrix.recDate = sessInfo(sessInd).recDate;
    goNogo_LFP_Matrix.odorPair1 = sessInfo(sessInd).pair1;
    goNogo_LFP_Matrix.odorPair2 = sessInfo(sessInd).pair2;

    % load corrected TIMESTAMPS
    syncFile = fullfile(savedir,'syncTime.mat');
    syncTime = load(syncFile);
    IO_ts = syncTime.Npx_timeStamps.IO.Ts_Sync;
    IO_ts = syncTime.Npx_timeStamps.IO.Ts_Sync;
    LFP_ts = syncTime.Npx_timeStamps.imec0.LFP_Ts_Sync;
    clear syncTime;

    %% start processing each probe
    for m = 1:length(sessInfo(sessInd).probeID)
        probeName =  sprintf('imec%d',sessInfo(sessInd).probeID(m));
        if p.savePlot
            savedir_plot = sprintf('%s%s%s%s%s%s',p.saveDir,'\Figures\',sessInfo(sessInd).animalID,'-date-',sessInfo(sessInd).recDate,'\LFP_Visualization');
            if ~exist(savedir_plot, 'dir')
                mkdir(savedir_plot);
            end
        end

        % load LFP
        % search for LFP file in the CAT folder, which is subsampled
        % filtered LFP
        folderTemp = dir(fullfile(sessInfo(sessInd).npx_path,strcat('catgt*_',probeName)));
        folderName = fullfile(folderTemp.folder,folderTemp.name);
        % detect LFP files
        lfpFile = dir(fullfile(folderName,'*.lf.bin'));
        metaFile = dir(fullfile(folderName,'*.lf.meta'));
        lfpFileName = fullfile(folderName,lfpFile.name);
        metaFileName = fullfile(folderName,metaFile.name);

        % Parse the corresponding metafile
        meta = ReadMeta(metaFileName);
        converter = ADconvert(meta); % convert to Volt
        % Get LFP
        lfp = [];
        for k = 1:length(LFP_ch)
            lfpTemp = ReadBinChannel(lfpFileName,meta,LFP_ch(k))*converter.AP(1)*10^6; % unit uV;
            lfp = [lfp;lfpTemp'];
        end

        %%
        % read in behavior time
        behaveFile = fullfile(savedir,'behaviorLabel.mat');
        load(behaveFile); % output: behaviorLabel

        if ~isempty(sessInfo(sessInd).pair1)

            trialNum = behaviorLabel.trialInfo_pair1.num;
            onsetInd = behaviorLabel.trialInfo_pair1.onsetInd;
            offsetInd = behaviorLabel.trialInfo_pair1.offsetInd;
            odor = behaviorLabel.trialInfo_pair1.odor; % 1,2,3,4
            odorNum = length(unique(odor));
            successLabel = behaviorLabel.trialInfo_pair1.success; % 0, premature lick, 1, successful trial
            outcomeLabel = behaviorLabel.trialInfo_pair1.outcome; % 1: Hit, 2: Miss, 3: FalseAlarm, 4: CorrectRejection

            trialInd_success = (successLabel==1);
            trialNum_success = sum(trialInd_success);
            trialNum_success = sum(trialInd_success);
            onsetInd_success = onsetInd(trialInd_success);
            offsetInd_success = offsetInd(trialInd_success);
            odorSuccess = odor(trialInd_success);
            outcomeSuccess = outcomeLabel(trialInd_success);

            LFP_delay = zeros(trialNum_success,ceil(p.delay_len*2500));

            for n = 1:trialNum_success
                h = figure;
                h.Position = [100,100,1600,1200];

                % plot LFP
                subplot(3,1,[1,2])
                plotInd = LFP_ts >= IO_ts(offsetInd(n)) & LFP_ts <= (IO_ts(offsetInd(n))+p.delay_len);
                for k = 1:length(LFP_ch)
                    plot(LFP_ts(plotInd),lfp(k,plotInd)-p.plotStep*k,'k','LineWidth',1)
                    hold on
                end
                axis tight
                % LFP_delay(n,:) = lfp(k,plotInd);

                TITLE1 = sprintf('%s-%s-pair1-onebox-IO-signal',behaviorLabel.animalID,behaviorLabel.recDate);
                TITLE2 = sprintf('Sucees-trial-%d-ordor-%d-outcome-%d',n,odorSuccess(n),outcomeSuccess(n));
                title({TITLE1;TITLE2},'Interpreter','None');

                % plot power spectrum
                subplot(3,1,3)
                wave.data = lfp(3,plotInd)';
                wave.timestamps = LFP_ts(plotInd);
                wave.samplingRate = 2500;
                spec = PowerSpectrum_Wavelet(wave,p.nCycles,p.dt,p.frange,'nfreqs',100);
                % freqInd = spec.freqs > 1 & spec.freqs < 100;
                imagesc(spec.timestamps,log2(spec.freqs),spec.amp')
                title('Spectrogram')
                SpecColorRange(spec.amp);
                colormap jet
                LogScale('y',2)
                ylabel('f (Hz)')
                xlabel('time (sec)')
                axis xy

                % subplot(4,1,4)
                % freqInd = spec.freqs >100;
                % imagesc(spec.timestamps,log2(spec.freqs(freqInd)),spec.amp(freqInd,:)')
                % % imagesc(spec.timestamps,log2(spec.freqs),spec.amp')
                % title('Spectrogram')
                % SpecColorRange(spec.amp(freqInd,:));
                % colormap jet
                % LogScale('y',2)
                % ylabel('f (Hz)')
                % xlabel('time (sec)')
                % axis xy
            

                figName = sprintf('%s%s%s-%s%s%d',savedir_plot,'\',behaviorLabel.animalID,behaviorLabel.recDate,'-LFP-Delay-successTrial-',n);
                print(figName,'-dpng','-r200');

                close all
            end

        end
    end

    if p.saveFile
        save(fullfile(savedir,'behaviorLabel.mat'), 'goNogo_LFP_Matrix','-v7.3');
    end
    fprintf('Finish behavioral analysis: %d\n', sessInd);
    clear goNogo_LFP_Matrix

end
close all
end
