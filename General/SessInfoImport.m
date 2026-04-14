function sessInfo = SessInfoImport(inFile)
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
read_total = sprintf('A5:W%d',max_number_of_i);
[~, ~, raw] = xlsread(inFile, worksheet ,read_total);
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};

% Transform the raw input into the right format
cellayerCH = cell(1,max(SessIdx));
diode = cell(1,max(SessIdx));

for row = 1:length(SessIdx_Raw)
    if ~isnan(SessIdx_Raw{row})  %only pick rows with valid i numbers
        i = SessIdx_Raw{row};
        
        if isnumeric(raw{row,2})
            ratID(i) = raw(row,2);
        end
        if ischar(raw{row,3});
            recDate{i} = raw{row,3};
        end
        
        
        if ischar(raw{row,4}); 
            neuropixel{i} = raw{row,4};
            oldFolder = cd;
%             try cd(neuropixel{i})  % check if paths are working
%             catch, errordlg(sprintf('Path incorrect! Cannot read i %d, path: %s', i, neuropixel{i})); end
%             cd(oldFolder);
        end
        
        if ischar(raw{row,5}); 
            IO{i} = raw{row,5};
            oldFolder = cd;
%             try cd(NIDQ{i})  % check if paths are working
%             catch, errordlg(sprintf('Path incorrect! Cannot read i %d, path: %s', i, NIDQ{i})); end
%             cd(oldFolder);
        end

        if isnumeric(raw{row,6})
            probeID{i} = {num2str(raw{row,6})};
        else
            probeIDTemp = strrep(raw{row,6},' ','');
            probeIDTemp = strrep(probeIDTemp,'''','');
            probeID{i} = strsplit(probeIDTemp,',');
        end

        if isnumeric(raw{row,7})
            NrxsynCh(i) = raw(row,7);
        end  
        if isnumeric(raw{row,8})
            NIsynCh(i) = raw(row,8);
        end  
        if isnumeric(raw{row,9})
            CamerasynCh(i) = raw(row,9);
        end 
        
        % phaseSess = strrep(raw{row,12},' ','');
        % phaseSess = strrep(phaseSess,'''','');
        % bSess{i} = strsplit(phaseSess,',');   
        % 
        % sleepSess = strrep(raw{row,13},' ','');
        % sleepSess = strrep(sleepSess,'''','');
        % sSess{i} = strsplit(sleepSess,','); 
        % 
        % openSess = strrep(raw{row,14},' ','');
        % openSess = strrep(openSess,'''','');
        % oSess{i} = strsplit(openSess,',');
        % 
        % EEGch =  strrep(raw{row,15},' ','');
        % EEGch = strrep(EEGch,'''','');
        % EEGDepth{i} = strsplit(EEGch,';');
        % 
        % Loc =  strrep(raw{row,16},' ','');
        % Loc = strrep(Loc,'''','');
        % EEGLoc{i} = strsplit(Loc,',');
        % 
        % EEGch2 =  strrep(raw{row,17},' ','');
        % EEGch2 = strrep(EEGch2,'''','');
        % EEGDepth2{i} = strsplit(EEGch2,';');
        % 
        % Loc2 =  strrep(raw{row,18},' ','');
        % Loc2 = strrep(Loc2,'''','');
        % EEGLoc2{i} = strsplit(Loc2,',');
        % 
        % probeid_hpc(i) = raw(row,19);
        % probeid_pfc(i) = raw(row,20);
        % probeid_mec(i) = raw(row,21);
    end
end

fld = {
    'session'
    'ratID'
    'Date'
    % 'neuralynx'
    'neuropixel'
    'IO'
    % 'diode'
    % 'NeuralynxSynCh'
    'probeID'
    'NrxSynCh'
    'NISynCh'
    'CameraSynCh'
    % 'ArdNISynCh'
    % 'sessDirs'
    % 'sleepDirs'
    % 'openDirs'
    % 'EEGDepth'
    % 'EEGLoc'
    % 'EEGDepth2'
    % 'EEGLoc2'
    % 'probeid_hpc'
    % 'probeid_pfc'
    % 'probeid_mec'
    };

vars = {
    num2cell(1:max(SessIdx))
    ratID
    recDate
    % neuralynx
    neuropixel
    IO
    % diode
    % synCh
    probeID
    NrxsynCh
    NIsynCh
    CamerasynCh
    % ArdNIsynCh
    % bSess
    % sSess
    % oSess
    % EEGDepth
    % EEGLoc
    % EEGDepth2
    % EEGLoc2
    % probeid_hpc
    % probeid_pfc
    % probeid_mec
    };

args = cell(length(fld)*2, 1);
args(1:2:end) = fld;
args(2:2:end) = vars;

sessInfo = struct(args{:});
sessInfo = sessInfo(:);

end