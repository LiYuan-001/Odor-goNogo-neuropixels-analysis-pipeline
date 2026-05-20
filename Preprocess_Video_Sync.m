% this code is to sync imec, obx (or NIDQ), and assign timestamps to each
% exsiting datapoint in ap.bin, lf.bin, and obx.bin (or NIDQ.bin)
% Li YUAN, 2026-Feb
function Preprocess_Video_Sync(inFile,analyzeSes)
close all
p.saveFile = 1;
p.saveDir = 'N:\yuan\AnalysisResults'; % change it to your own folder

p.dataType = 'D';
% Read in input information
sessInfo = SessInfoImport_Npx(inFile);

for sessInd = analyzeSes(1:end)

    close all
    
    videoTime.animalID = sessInfo(sessInd).animalID;
    videoTime.recDate = sessInfo(sessInd).recDate;

    p.videoFs = sessInfo(sessInd).video_rate;

    if p.saveFile
        savedir = sessInfo(sessInd).npx_path;
        % savedir = sprintf('%s%s%s%s%s',p.saveDir,'\ProcessedData\',sessInfo(sessInd).animalID,'-date-',sessInfo(sessInd).recDate); 
        % if ~exist(savedir, 'dir')
        %     mkdir(savedir);
        % end
    end
    p.cam_frame_ch = sessInfo(sessInd).cam_frame_ch;

    % load(fullfile(savedir,'syncTime.mat'));
    % IO_type = Npx_timeStamps.IO.Type;
    % 
    % % Check if either file exists
    % if strcmpi(IO_type,'nidq')
    %     IOFile = dir(fullfile(sessInfo(sessInd).npx_path,'*.nidq.bin'));
    %     metaFile = dir(fullfile(sessInfo(sessInd).npx_path,'*.nidq.meta'));
    % elseif strcmpi(IO_type,'obx')
    %     IOFile = dir(fullfile(sessInfo(sessInd).npx_path,'*.obx.bin'));
    %     metaFile = dir(fullfile(sessInfo(sessInd).npx_path,'*.obx.meta'));        
    % else
    %     fprintf('Error,No NIDQ.meta or obx.meta found in %s\n', folderPath);
    %     return
    % end
    % 
    % IOFileName = IOFile.name;
    % metaFileName = metaFile.name;

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



    % Parse the corresponding metafile
    meta_IO = ReadMeta(fullfile(sessInfo(sessInd).npx_path,metaFileName));
    IO_Fs = SampRate(meta_IO);
    % Get first one channel of LFP
    dataArray = ReadBin(fullfile(sessInfo(sessInd).npx_path,IOFileName),meta_IO);

    if contains(IO_type,'nidq')
        if p.dataType == 'A'
            dataArray = GainCorrectNI(dataArray, [0:7], meta_IO); % change the [0:7] is only part of channels activated in setting
            camTTL = dataArray(p.cam_frame_ch,:);
            camTTL = camTTL > 2;
        else
            % (1-based). For imec data there is never more than one saved digital word.
            dw = 1;
            % Read these lines in dw (0-based).
            dLineList = [0:7]; % The 2090A NIDQ
            digArray = ExtractDigital(dataArray, meta_IO, dw, dLineList);
            digArray2 = double(digArray);
            clear digArray
            camTTL = digArray2(p.cam_frame_ch,:); % onebox specific
            clear digArray2
        end

    elseif contains(IO_type,'obx')
        if p.dataType == 'A'
            dataArray = GainCorrectNI(dataArray, [0:11], meta_IO); % change the [0:11] is only part of channels activated in setting
            camTTL = dataArray(p.cam_frame_ch,:);
            camTTL = camTTL > 2;
        else
            % (1-based). For imec data there is never more than one saved digital word.
            dw = 1;
            % Read these lines in dw (0-based).
            dLineList = [0:11];
            digArray = ExtractDigital(dataArray, meta_IO, dw, dLineList);
            digArray2 = double(digArray);
            clear digArray
            camTTL = digArray2(p.cam_frame_ch,:); % onebox specific
            clear digArray2
        end
    end

    clear dataArray

    camTTL = camTTL(:);
    camTTL_diff = [0;diff(camTTL)];
    droppingInd = camTTL_diff==-1;
    framidx = find(droppingInd);
    camTTL_frameCount = sum(droppingInd);

    %% read video
    p.imageFolder = sessInfo(sessInd).video_path;
    % Make a video from images in p.imageFolder using image-header timestamps
    % Assumes each image file contains a readable timestamp in its metadata.
    % Output video frame rate is p.videoFs.

    imgExt = '*.tif';   % change if needed, e.g. '*.png' or '*.jpg'
    files = dir(fullfile(p.imageFolder, imgExt));
    if isempty(files)
        error('No images found in %s', p.imageFolder);
    end

    % Sort by filename
    [~, idx] = sort({files.name});
    files = files(idx);

    nFile = numel(files);
    imgTime = NaT(nFile,1);

    % compare numbers of frames
    if camTTL_frameCount == nFile
        fprintf('Frames match TTL\n')
        videoTime.idx = framidx;
        % videoTime.ts = Npx_timeStamps.IO.Ts_Sync(framidx);
    % elseif camTTL_frameCount == nFile+1 % this comes from last frame is triggered but not saved
    %     fprintf('Frames match TTL\n')
    %     videoTime.ind = framidx(1:end-1);
    %     videoTime.ts = Npx_timeStamps.IO.Ts_Sync(framidx(1:end-1));
    else

        fprintf('Frames doesnt match TTL\n')
        break
    end


    % Output file
    outFile = fullfile(sessInfo(sessInd).npx_path, 'output_video.avi');
    % delete existing file, otherwise concatenate
    if exist(outFile,'file')
        delete(outFile)
        fprintf('Deleted exisiting video file.\n')
    end
    
    v = VideoWriter(outFile, 'Motion JPEG AVI');
    v.FrameRate = p.videoFs;
    open(v);

    nextProgress = 10;

    % Read timestamps from image headers
    for k = 1:nFile
        fname = fullfile(files(k).folder, files(k).name);
        info = imfinfo(fname);

        % Try common timestamp fields
        if isfield(info, 'DateTime') && ~isempty(info.DateTime)
            imgTime(k) = datetime(info.DateTime, 'InputFormat', 'yyyy:MM:dd HH:mm:ss');
        elseif isfield(info, 'ImageDescription') && ~isempty(info.ImageDescription)
            % Example: parse a timestamp embedded in ImageDescription
            txt = string(info.ImageDescription);

            % Edit this regexp to match your actual timestamp format
            token = regexp(txt, '\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}(\.\d+)?', 'match', 'once');
            if ~isempty(token)
                imgTime(k) = datetime(token, 'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSS');
            else
                error('No recognizable timestamp in ImageDescription for %s', files(k).name);
            end
        else
            error('No supported timestamp field found in %s', files(k).name);
        end

        % read image and write video
        img = imread(fname);
        writeVideo(v, img);
        % Show progress every 10%
        progress = floor(k / nFile * 100);

        if progress >= nextProgress
            fprintf('Progress: %d%%\n', nextProgress);
            nextProgress = nextProgress + 10;
        end
    end

    fprintf('Progress: 100%% Done.\n');
    close(v);
    fprintf('Saved video to: %s\n', outFile);

    % Convert to seconds relative to first frame
    tSec = seconds(imgTime - imgTime(1));
    % Target uniform video times
    % tVideo = (0 : 1/p.videoFs : (length(tSec)-1)/p.videoFs)';
    videoTime.tVideo_Origin = tSec;

    if p.saveFile
        save(fullfile(savedir,'videoTime.mat'), 'videoTime','-v7.3');
    end

    clear videoTime
    fprintf('Finish Imec video sync\n');

end
close all
end
