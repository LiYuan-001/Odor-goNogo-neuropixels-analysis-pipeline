function sessInfo = SessInfoImport_Nrpix(inFile)
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
            animalID(i) = raw(row,2);
        end
        if ischar(raw{row,3})
            recDate{i} = raw{row,3};
        end
              
        if ischar(raw{row,4})
            nrpix_path{i} = raw{row,4};
            % oldFolder = cd;
            % try cd(neuropixel{i})  % check if paths are working
            % catch, errordlg(sprintf('Path incorrect! Cannot read i %d, path: %s', i, neuropixel{i})); end
            % cd(oldFolder);
        end

        if ischar(raw{row,5})
            nrpix_IO_path{i} = raw{row,5};
            % oldFolder = cd;
            % try cd(neuropixel{i})  % check if paths are working
            % catch, errordlg(sprintf('Path incorrect! Cannot read i %d, path: %s', i, neuropixel{i})); end
            % cd(oldFolder);
        end
        
        if ischar(raw{row,6})
            NI_path{i} = raw{row,6};
            % oldFolder = cd;
            % try cd(NIDQ{i})  % check if paths are working
            % catch, errordlg(sprintf('Path incorrect! Cannot read i %d, path: %s', i, NIDQ{i})); end
            % cd(oldFolder);
        end
             
        if isnumeric(raw{row,7})
            nrpix_synCh(i) = raw(row,7);
        end  
        if isnumeric(raw{row,8})
            nrpix_IO_synCh(i) = raw(row,8);
        end  
        if isnumeric(raw{row,9})
            NI_synCh(i) = raw(row,9);
        end  
        if isnumeric(raw{row,10})
            ArdNI_synCh(i) = raw(row,10);
        end 
        
        odorSess = strrep(raw{row,11},' ','');
        odorSess = strrep(odorSess,'''','');
        oSess{i} = strsplit(odorSess,',');   
        
        sleepSess = strrep(raw{row,12},' ','');
        sleepSess = strrep(sleepSess,'''','');
        sSess{i} = strsplit(sleepSess,','); 

        EEGch =  strrep(raw{row,13},' ','');
        EEGch = strrep(EEGch,'''','');
        EEGDepth{i} = strsplit(EEGch,';');
        
        Loc =  strrep(raw{row,14},' ','');
        Loc = strrep(Loc,'''','');
        EEGLoc{i} = strsplit(Loc,',');
        
        EEGch2 =  strrep(raw{row,15},' ','');
        EEGch2 = strrep(EEGch2,'''','');
        EEGDepth2{i} = strsplit(EEGch2,';');
        
        Loc2 =  strrep(raw{row,16},' ','');
        Loc2 = strrep(Loc2,'''','');
        EEGLoc2{i} = strsplit(Loc2,',');

        probeID = strrep(raw{row,17},' ','');
        probeID = strrep(probeID,'''','');
        probeID_sync{i} = strsplit(probeID,',');
        
        probeID_hpc(i) = raw(row,18);
        probeID_pfc(i) = raw(row,19);
    end
end

fld = {
    'session'
    'animalID'
    'Date'
    'Nrpix_path'
    'Nrpix_IO_path'
    'NIDQ_path'
    'Nrpix_IO_SynCh'
    'Nrpix_SynCh'
    'Nrpix_IO_SynCh'
    'Ard_SynCh'
    'Odor'
    'sleepDirs'
    'EEGDepth'
    'EEGLoc'
    'EEGDepth2'
    'EEGLoc2'
    'probeID_Syn'
    'probeID_hpc'
    'probeID_pfc'
    };

vars = {
    num2cell(1:max(SessIdx))
    animalID
    recDate
    nrpix_path
    nrpix_IO_path
    NI_path
    nrpix_synCh
    nrpix_IO_synCh
    NI_synCh
    ArdNI_synCh
    oSess
    sSess
    EEGDepth
    EEGLoc
    EEGDepth2
    EEGLoc2
    probeID_sync
    probeID_hpc
    probeID_pfc
    };

args = cell(length(fld)*2, 1);
args(1:2:end) = fld;
args(2:2:end) = vars;

sessInfo = struct(args{:});
sessInfo = sessInfo(:);

end