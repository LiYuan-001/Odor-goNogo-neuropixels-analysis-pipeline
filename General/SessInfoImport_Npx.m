function sessInfo = SessInfoImport_Npx(inFile)
% Define input parameters
% number of rows the program reads, increase if more then 700 rows are used up
max_number_of_i = 5000;
worksheet = 'Sheet1';

% Read first column, to define all i numbers that occur
read_i = sprintf('A5:A%d',max_number_of_i);
[~, ~, SessIdx_Raw] = xlsread(inFile, worksheet, read_i);
SessIdx = reshape([SessIdx_Raw{:,1}],size(SessIdx_Raw));
SessIdx(isnan(SessIdx)) = [];

% Now import the rest of the data
read_total = sprintf('A5:AZ%d',max_number_of_i);
[~, ~, raw] = xlsread(inFile, worksheet ,read_total);
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};

% Transform the raw input into the right format

for row = 1:length(SessIdx)
    if ~isnan(SessIdx_Raw{row})  %only pick rows with valid i numbers
        sessInd = SessIdx_Raw{row};

        if isnumeric(raw{row,2})
            animalID{sessInd} = char(string(raw{row,2}));
        else
            animalID{sessInd} = char(raw{row,2});
        end

        if ~isempty(raw{row,3})
            recDate{sessInd} = string(raw{row,3});
        else
            recDate{sessInd} = "";
        end

        if ischar(raw{row,4})
            npx_path{sessInd} = raw{row,4};
            % oldFolder = cd;
            % try cd(npx_path{sessInd})  % check if paths are working
            % catch, errordlg(sprintf('Path incorrect! Cannot read i %d, path: %s', i, npx_path{sessInd})); end
            % cd(oldFolder);
        end

        if ischar(raw{row,5})
            video_path{sessInd} = raw{row,5};
            % oldFolder = cd;
            % try cd(video_path{sessInd})  % check if paths are working
            % catch, errordlg(sprintf('Path incorrect! Cannot read i %d, path: %s', i, video_path{sessInd})); end
            % cd(oldFolder);
        end

        video_rate(sessInd) = raw(row,6);

        probeID_raw = raw{row,7};
        if isnumeric(probeID_raw)
            probeID_temp = probeID_raw;
        elseif ischar(probeID_raw) || isstring(probeID_raw)
            txt = char(string(probeID_raw));
            txt = strrep(txt, ' ', '');   % remove spaces
            parts = regexp(txt, '[,;]', 'split');
            probeID_temp = str2double(parts);
            probeID_temp = probeID_temp(~isnan(probeID_temp));
        else
            probeID_temp = [];
        end
        probeID{sessInd} = probeID_temp;

        IO_synCh(sessInd) = raw(row,8);

        % correct from machine to matlab 
        % matlab count from 1. machine count from 0

        if isnumeric(raw{row,9})
            splitIndicator(sessInd) =  num2cell(raw{row,9} + 1);
        else
            splitIndicator(sessInd) =  raw(row,9);
        end

        % splitTemp =  strrep(raw{row,10},' ','');
        % splitTemp = strrep(splitTemp,'''','');
        % splitTemp = lower(splitTemp);
        % splitNames{sessInd} = regexp(splitTemp, '[,;]', 'split');

        if isnumeric(raw{row,10})
            cam_frame_ch(sessInd)= num2cell(raw{row,10} + 1);
        else
            cam_frame_ch(sessInd)= raw(row,10);
        end

        if isnumeric(raw{row,11})
            odorA_ch(sessInd)= num2cell(raw{row,11} + 1);
        else
            odorA_ch(sessInd)= raw(row,11);
        end

        if isnumeric(raw{row,12})
            odorB_ch(sessInd)= num2cell(raw{row,12} + 1);
        else
            odorB_ch(sessInd)= raw(row,12);
        end

        if isnumeric(raw{row,13})
            odorC_ch(sessInd)= num2cell(raw{row,13} + 1);
        else
            odorC_ch(sessInd)= raw(row,13);
        end

        if isnumeric(raw{row,14})
            odorD_ch(sessInd)= num2cell(raw{row,14} + 1);
        else
            odorD_ch(sessInd)= raw(row,14);
        end

        if isnumeric(raw{row,15})
            sucrose_ch(sessInd)= num2cell(raw{row,15} + 1);
        else
            sucrose_ch(sessInd)= raw(row,15);
        end

        if isnumeric(raw{row,16})
            quinine_ch(sessInd)= num2cell(raw{row,16} + 1);
        else
            quinine_ch(sessInd)= raw(row,16);
        end

        if isnumeric(raw{row,17})
            led_ch(sessInd)= num2cell(raw{row,17} + 1);
        else
            led_ch(sessInd)= raw(row,17);
        end

        if isnumeric(raw{row,18})
            lick_ch(sessInd)= num2cell(raw{row,18} + 1);
        else
            lick_ch(sessInd)= raw(row,18);
        end

        if isnumeric(raw{row,19})
            stim_ch(sessInd)= num2cell(raw{row,19} + 1);
        else
            stim_ch(sessInd)= raw(row,19);
        end

        odor_raw = raw{row,20};
        pair1{sessInd} = num2cell((odor_raw));

        odor_raw = raw{row,21};
        pair2{sessInd} = num2cell((odor_raw));

        delay_len(sessInd) = raw(row,22);
        response_len(sessInd) = raw(row,23);
        ITI(sessInd) = raw(row,24);
        ihiTag = raw(row,25);

        sleepSess = strrep(raw{row,26},' ','');
        sleepSess = strrep(sleepSess,'''','');
        sblock{sessInd} = regexp(sleepSess, '[,;]', 'split');

        EEGch =  strrep(raw{row,27},' ','');
        EEGch = strrep(EEGch,'''','');
        EEGDepth{sessInd} = regexp(EEGch, '[,;]', 'split');

        Loc =  strrep(raw{row,28},' ','');
        Loc = strrep(Loc,'''','');
        Loc = lower(Loc);
        EEGLoc{sessInd} = regexp(Loc, '[,;]', 'split');

        EEGch2 =  strrep(raw{row,29},' ','');
        EEGch2 = strrep(EEGch2,'''','');
        EEGDepth2{sessInd} = regexp(EEGch2, '[,;]', 'split');

        Loc2 =  strrep(raw{row,30},' ','');
        Loc2 = strrep(Loc2,'''','');
        Loc2 = lower(Loc2);
        EEGLoc2{sessInd} = regexp(Loc2, '[,;]', 'split');

        % probeID = strrep(raw{row,17},' ','');
        % probeID = strrep(probeID,'''','');
        % probeID_sync{i} = strsplit(probeID,',');
        %
        % probeID_hpc(i) = raw(row,18);
        % probeID_pfc(i) = raw(row,19);
    end
end

fld = {
    'sessID'
    'animalID'
    'recDate'
    'npx_path'
    'video_path'
    'video_rate'
    'probeID'
    'IO_synCh'
    'splitIndicator'
    'cam_frame_ch'
    'odorA_ch'
    'odorB_ch'
    'odorC_ch'
    'odorD_ch'
    'sucrose_ch'
    'quinine_ch'
    'led_ch'
    'lick_ch'
    'stim_ch'
    'pair1'
    'pair2'
    'delay_len'
    'response_len'
    'ITI'
    'ihiTag'
    'sblock'
    'EEGDepth'
    'EEGLoc'
    'EEGDepth2'
    'EEGLoc2'
    % 'probeID_Syn'
    % 'probeID_hpc'
    % 'probeID_pfc'
    };

vars = {
    num2cell(1:max(SessIdx))
    animalID
    recDate
    npx_path
    video_path
    video_rate
    probeID
    IO_synCh
    splitIndicator
    cam_frame_ch
    odorA_ch
    odorB_ch
    odorC_ch
    odorD_ch
    sucrose_ch
    quinine_ch
    led_ch
    lick_ch
    stim_ch
    pair1
    pair2
    delay_len
    response_len
    ITI
    ihiTag
    sblock
    EEGDepth
    EEGLoc
    EEGDepth2
    EEGLoc2
    };

args = cell(length(fld)*2, 1);
args(1:2:end) = fld;
args(2:2:end) = vars;

sessInfo = struct(args{:});
sessInfo = sessInfo(:);

end