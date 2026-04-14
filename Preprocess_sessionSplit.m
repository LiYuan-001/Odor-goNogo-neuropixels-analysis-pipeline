% this code is to get the index where the recording big session split into
% small blocks, such as pair1, pair2, sleep etc
% Li YUAN, 2026-Mar
function Preprocess_sessionSplit(inFile,analyzeSes)
close all
p.saveFile = 1;
p.saveDir = 'N:\yuan\AnalysisResults'; % change it to your own folder
p.dataType = 'A';
% Read in input information
sessInfo = SessInfoImport_Npx(inFile);

for sessInd = analyzeSes(1:end)

    close all

    if p.saveFile
        savedir = sprintf('%s%s%s%s%s',p.saveDir,'\ProcessedData\',sessInfo(sessInd).animalID,'-date-',sessInfo(sessInd).recDate); 
        % if ~exist(savedir, 'dir')
        %     mkdir(savedir);
        % end
    end

    behaviorBlocks.animalID = sessInfo(sessInd).animalID;
    behaviorBlocks.recDate = sessInfo(sessInd).recDate;
    behaviorBlocks.odorPair1 = sessInfo(sessInd).pair1;
    behaviorBlocks.odorPair2 = sessInfo(sessInd).pair2;

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
        event.npx.file = 1;
        % Parse the corresponding metafile
        meta_IO = ReadMeta(fullfile(sessInfo(sessInd).npx_path,metaFileName));
        IO_Fs = SampRate(meta_IO);
        % Get first one channel of LFP
        dataArray = ReadBin(fullfile(sessInfo(sessInd).npx_path,IOFileName),meta_IO);

        if contains(IO_type,'nidq')
            if p.dataType == 'A'
                dataArray = GainCorrectNI(dataArray, [0:7]+1, meta_IO); % change the [0:7] is only part of channels activated in setting
                behaveData = dataArray;
            else
                % (1-based). For imec data there is never more than one saved digital word.
                dw = 1;
                % Read these lines in dw (0-based).
                dLineList = [0:7] + 1; % correct counts for matlab
                digArray = ExtractDigital(dataArray, meta_IO, dw, dLineList);
                digArray2 = double(digArray);
                clear digArray
                behaveData = digArray2;
                clear digArray2
            end

        elseif contains(IO_type,'obx')
            if p.dataType == 'A'
                dataArray = GainCorrectOBX(dataArray, [0:11]+1, meta_IO); % change the [0:11] is only part of channels activated in setting
                behaveData = dataArray;
            else
                % (1-based). For imec data there is never more than one saved digital word.
                dw = 1;
                % Read these lines in dw (0-based).
                dLineList = [0:11] + 1; % correct counts for matlab
                digArray = ExtractDigital(dataArray, meta_IO, dw, dLineList);
                digArray2 = double(digArray);
                clear digArray
                behaveData = digArray2;
                clear digArray2
            end
        end
        clear dataArray

        IO_ts = (0:size(behaveData,2)-1) / IO_Fs;

        splitIndicator = sessInfo(sessInd).splitIndicator;
        % splitNames = sessInfo(sessInd).splitNames;
        % splitInd = find(splitNames == 'odor');
        switchSig = behaveData(splitIndicator,:) > 4;
        changeInd = find(diff(switchSig) ==1);

        eventTable = annotateSwitchEvents(switchSig, IO_Fs);
        eventStruct = table2struct(eventTable,'ToScalar', true);
        behaviorBlocks.eventNames = eventStruct.eventName;
        behaviorBlocks.signalInd = eventStruct.signalIdx;

    if p.saveFile
        save(fullfile(savedir,'behaviorBlocks.mat'), 'behaviorBlocks','-v7.3');
    end
    fprintf('Finish behavioral block analysis: %d\n', sessInd);
    clear behaviorBlocks

end
close all
end
