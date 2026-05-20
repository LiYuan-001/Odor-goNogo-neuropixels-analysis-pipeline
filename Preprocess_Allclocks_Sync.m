% this code is to sync imec, obx (or NIDQ), and assign timestamps to each
% exsiting datapoint in ap.bin, lf.bin, and obx.bin (or NIDQ.bin)
% Li YUAN, 2026-Feb
function Preprocess_Allclocks_Sync(inFile,analyzeSes)
close all
p.saveFile = 1;
p.saveDir = 'N:\yuan\AnalysisResults'; % change it to your own folder

% Read in input information
sessInfo = SessInfoImport_Npx(inFile);

for sessInd = analyzeSes(1:end)

    close all

    sessionName = sprintf('%s%s',sessInfo(sessInd).animalID,sessInfo(sessInd).recDate);
    % get sync channel ID
    p.npx_synCh = 385;
    p.IO_synCh = sessInfo(sessInd).IO_synCh;

    if p.saveFile
        savedir = sprintf('%s%s%s%s%s',p.saveDir,'\ProcessedData\',sessInfo(sessInd).animalID,'-date-',sessInfo(sessInd).recDate); 
        if ~exist(savedir, 'dir')
            mkdir(savedir);
        end
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
        fprintf('Error,No IO.meta or obx.meta found in %s\n', folderPath);
        return
    end

    % Parse the corresponding metafile
    meta_IO = ReadMeta(fullfile(sessInfo(sessInd).npx_path,metaFileName));
    IO_Fs = SampRate(meta_IO);
    % % read in data 
    % dataArray = ReadBin(fullfile(sessInfo(sessInd).npx_path,IOFileName),meta_IO);
    % 
    % if contains(IO_type,'nidq')
    %     % (1-based). For imec data there is never more than one saved digital word.
    %     dw = 1;
    %     % Read these lines in dw (0-based).
    %     dLineList = [0:7]; % The 2090A NIDQ
    %     digArray = ExtractDigital(dataArray, meta_IO, dw, dLineList);
    %     digArray2 = double(digArray);
    %     clear digArray
    %     IO_sync_1 = digArray2(sessInfo(sessInd).IO_synCh,:); % onebox specific
    %     clear digArray2
    % 
    % elseif contains(IO_type,'obx')
    %     IO_sync_1 = dataArray(sessInfo(sessInd).IO_synCh,:); % onebox specific
    % else
    %     error('No IO signal for sync')
    % end
    % clear dataArray

    IO_sync_1 = ReadBinChannel(fullfile(sessInfo(sessInd).npx_path,IOFileName),meta_IO,sessInfo(sessInd).IO_synCh);

    % Calculate and assign time
    IO_Ts = (1:length(IO_sync_1))./IO_Fs;
    % use IO as main clock
    Npx_timeStamps.IO.Type = IO_type;
    Npx_timeStamps.IO.Fs = IO_Fs;
    Npx_timeStamps.IO.Ts_Sync = IO_Ts;

    
    figure(1)
    plot(IO_Ts,IO_sync_1,'k')
    text(0,0,'IO channel','Interpreter','None');
    hold on

    fprintf('Finish IO file. Waiting for probe timestamps Process\n');

    %% Neuropixels
    for m = 1:length(sessInfo(sessInd).probeID) % in case multiple probes
        probeName = sprintf('imec%d',sessInfo(sessInd).probeID(m));
        probeName_subFolder = sprintf('*_imec%d',sessInfo(sessInd).probeID(m));
        % detect file name
        %     oldFolder = cd;
        %     cd(sessInfo(i).neuropixel);
        folderTemp = dir(fullfile(sessInfo(sessInd).npx_path,probeName_subFolder));
        folderName = [];
        for k = 1:length(folderTemp)
            folderName{k} = fullfile(folderTemp(k).folder,folderTemp(k).name);
        end
        %     cd(folderName(1).name);

        % go to folder with both LF and AP file
        % detect file name
        % since neuropixel 2.0, LFP file is not generated in spikeGLX but in
        % the catGT. I edit code here to flexibablly detects LFP files in any
        % folder  ----- Li Yuan, Jan-2024
        lfpFile = [];
        metaFile =[];
        lfFileInd = [];
        apFile =[];
        ap_metaFile =[];
        apFileInd = [];
        for k = 1:length(folderTemp)
            lfpFile{k} = dir(fullfile(folderName{k},'*.lf.bin'));
            metaFile{k} = dir(fullfile(folderName{k},'*.lf.meta'));
            lfFileInd = [lfFileInd,~isempty(lfpFile{k})];

            apFile{k} = dir(fullfile(folderName{k},'*.ap.bin'));
            ap_metaFile{k} = dir(fullfile(folderName{k},'*.ap.meta'));
            apFileInd = [apFileInd,~isempty(apFile{k})];
        end

        lfpFileInd = find(lfFileInd==1);
        apFileInd = find(apFileInd==1);

        lfpFileName = fullfile(folderName{lfpFileInd(1)},lfpFile{lfpFileInd(1)}.name);
        metaFileName = fullfile(folderName{lfpFileInd(1)},metaFile{lfpFileInd(1)}.name);
        apFileName = fullfile(folderName{apFileInd(1)},apFile{apFileInd(1)}.name);
        ap_metaFileName = fullfile(folderName{apFileInd(1)},ap_metaFile{apFileInd(1)}.name);
        % route back to original folder
        %     cd(oldFolder);


        % Parse the corresponding metafile
        meta = ReadMeta(metaFileName);
        meta_ap = ReadMeta(ap_metaFileName);

        % get information for probe type
        % list of probe types with NP 1.0 imro format
        np1_imro = [0,1020,1030,1200,1100,1120,1121,1122,1123,1300];
        if isfield(meta_ap,'imDatPrb_type')
            probeType = str2double(meta_ap.imDatPrb_type);
        else
            probeType = 0;
        end
        Npx_timeStamps.(probeName).probeType = probeType;

        Npx_timeStamps.(probeName).LFP_Fs = SampRate(meta);
        Npx_timeStamps.(probeName).AP_Fs = SampRate(meta_ap);
        converter = ADconvert(meta_ap); % convert to Volt
        Npx_timeStamps.(probeName).converter = converter;


        % show meta in the command window
        % meta

        % Get original timestamps through index before syncing
        % LFP
        nChan = str2double(meta.nSavedChans);
        nFileSamp = str2double(meta.fileSizeBytes) / (2 * nChan);
        LFP_Ts_Origin = (0:nFileSamp-1)/Npx_timeStamps.(probeName).LFP_Fs;
        % AP
        nChan_ap = str2double(meta_ap.nSavedChans);
        nFileSamp_ap = str2double(meta_ap.fileSizeBytes) / (2 * nChan_ap);
        AP_Ts_Origin = (0:nFileSamp_ap-1)/Npx_timeStamps.(probeName).AP_Fs;

        % Nrp.(probeName).LFP_Ts_Sync = Nrp.(probeName).LFP_Ts_Origin;
        % Nrp.(probeName).AP_Ts_Sync = Nrp.(probeName).AP_Ts_Origin;

        % do sync
        lfp_sync = ReadBinChannel(lfpFileName,meta,p.npx_synCh);

        % detect rising egde in the sync_1 and align the IMEC time to IO time
        % Sometimes IMEC could have missing samples thats why using NIDA as
        % time
        % find exception (usually only 0)
        while any(lfp_sync>80)
            indTemp = find(lfp_sync>80);
            % assign exception to normal value (barely happen)
            lfp_sync(indTemp) = lfp_sync(indTemp-1);
        end
        % assign exception to normal value (barely happen)
        sync_Temp = (lfp_sync>50);
        Npx_timeStamps.(probeName).LFP_Ts_Sync = syncTime2(IO_sync_1,sync_Temp,Npx_timeStamps.IO.Ts_Sync,Npx_timeStamps.(probeName).LFP_Fs,1);
        % ap time stamp sync
        ap_sync = ReadBinChannel(apFileName,meta_ap,p.npx_synCh);

        % find exception
        %     ap_sync = detrend(ap_sync);
        while any(ap_sync>80)
            indTemp = find(ap_sync>80);
            % assign exception to normal value (barely happen)
            ap_sync(indTemp) = ap_sync(indTemp-5);
        end
        sync_Temp = (ap_sync>50);
        clear ap_sync
        Npx_timeStamps.(probeName).AP_Ts_Sync = syncTime2(IO_sync_1,sync_Temp,Npx_timeStamps.IO.Ts_Sync,Npx_timeStamps.(probeName).AP_Fs,1);

        plot(Npx_timeStamps.(probeName).LFP_Ts_Sync,lfp_sync + 200,'k')
        text(max(Npx_timeStamps.(probeName).LFP_Ts_Sync)-3,200,'Probe Ts sync','Interpreter','None');
        text(0,200,'Probe 0 Ts sync','Interpreter','None');
        hold on

        plot(LFP_Ts_Origin,lfp_sync + 100)
        text(max(Npx_timeStamps.(probeName).LFP_Ts_Sync)-3,100,'Probe Ts origin','Interpreter','None');
        text(max(Npx_timeStamps.(probeName).LFP_Ts_Sync)-3,0,'IO channel','Interpreter','None');
        text(0,100,'Probe 0 Ts origin','Interpreter','None');
        xlabel('Time (s)')
        xlim([max(Npx_timeStamps.(probeName).LFP_Ts_Sync)-3 max(Npx_timeStamps.(probeName).LFP_Ts_Sync)])
    end


    if p.saveFile
        save(fullfile(savedir,'syncTime.mat'), 'sessionName','Npx_timeStamps','-v7.3');
        save(fullfile(sessInfo(sessInd).npx_path,'syncTime.mat'), 'sessionName','Npx_timeStamps','-v7.3');
    end

    figName = sprintf('%s%s%s%s%s',savedir,'\Sync_illustration_',sessInfo(sessInd).animalID,'-Day-',sessInfo(sessInd).recDate);
    print(figName,'-dpng','-r300');

    clear Npx_timeStamps
    fprintf('Finish Imec sync for session %d\n', sessInd);

end
close all
end
