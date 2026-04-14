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
        savedir = sprintf('%s%s%s%s%s',p.saveDir,'\ProcessedData\',sessInfo(sessInd).animalID,'-date-',sessInfo(sessInd).recDate); 
        % if ~exist(savedir, 'dir')
        %     mkdir(savedir);
        % end
    end
    p.cam_frame_ch = sessInfo(sessInd).cam_frame_ch + 1; % matlab count from 1

    load(fullfile(savedir,'syncTime.mat'));
    IO_type = Npx_timeStamps.IO.Type;

    % Check if either file exists
    if strcmpi(IO_type,'nidq')
        IOFile = dir(fullfile(sessInfo(sessInd).npx_path,'*.nidq.bin'));
        metaFile = dir(fullfile(sessInfo(sessInd).npx_path,'*.nidq.meta'));
    elseif strcmpi(IO_type,'obx')
        IOFile = dir(fullfile(sessInfo(sessInd).npx_path,'*.obx.bin'));
        metaFile = dir(fullfile(sessInfo(sessInd).npx_path,'*.obx.meta'));        
    else
        fprintf('Error,No NIDQ.meta or obx.meta found in %s\n', folderPath);
        return
    end

    IOFileName = IOFile.name;
    metaFileName = metaFile.name;

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
        videoTime.ts = Npx_timeStamps.IO.Ts_Sync(droppingInd);
    else
        fprintf('Frames doesnt match TTL/n')
        break
    end

    % % Read timestamps from image headers
    % for k = 1:nFile
    %     fname = fullfile(files(k).folder, files(k).name);
    %     info = imfinfo(fname);
    % 
    %     % Try common timestamp fields
    %     if isfield(info, 'DateTime') && ~isempty(info.DateTime)
    %         imgTime(k) = datetime(info.DateTime, 'InputFormat', 'yyyy:MM:dd HH:mm:ss');
    %     elseif isfield(info, 'ImageDescription') && ~isempty(info.ImageDescription)
    %         % Example: parse a timestamp embedded in ImageDescription
    %         txt = string(info.ImageDescription);
    % 
    %         % Edit this regexp to match your actual timestamp format
    %         token = regexp(txt, '\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}(\.\d+)?', 'match', 'once');
    %         if ~isempty(token)
    %             imgTime(k) = datetime(token, 'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSS');
    %         else
    %             error('No recognizable timestamp in ImageDescription for %s', files(k).name);
    %         end
    %     else
    %         error('No supported timestamp field found in %s', files(k).name);
    %     end
    % end
    % 
    % % Convert to seconds relative to first frame
    % tSec = seconds(imgTime - imgTime(1));
    % % Target uniform video times
    % tVideo = (0 : 1/p.videoFs : (length(tSec)-1)/p.videoFs)';
    % videoTime.tVideo = tVideo;

    imgIdx = 1:nFile;
    % Output file
    outFile = fullfile(sessInfo(sessInd).npx_path, 'output_video.mp4');
    v = VideoWriter(outFile, 'MPEG-4');
    v.FrameRate = p.videoFs;
    open(v);

    for k = 1:numel(imgIdx)
        img = imread(fullfile(files(imgIdx(k)).folder, files(imgIdx(k)).name));

        % Convert grayscale to RGB for safer video writing
        if ndims(img) == 2
            img = repmat(img, 1, 1, 3);
        end

        writeVideo(v, img);
    end

    close(v);
    fprintf('Saved video to: %s\n', outFile);


    if p.saveFile
        save(fullfile(savedir,'videoTime.mat'), 'videoTime','-v7.3');
    end

    clear videoTime
    fprintf('Finish Imec video sync\n');

end
close all
end
