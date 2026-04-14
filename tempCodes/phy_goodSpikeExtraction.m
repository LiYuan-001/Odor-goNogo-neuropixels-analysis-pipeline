mouseID = 0001;
recDate = rec_20260101;
p.savePlot = 1;
if p.savePlot
    % directory for plot figures
    % generate a folder for each rat eah day under the current folder
    savedir_1 = sprintf('%s%s%d%s%s%s',cd,'\Figures\',mouseID,'-day-',recDate,'\SpikeShape\Imec0');
    if ~exist(savedir_1, 'dir')
        mkdir(savedir_1);
    end
    % delete(strcat(savedir_1,'\*'));
end

p.writeToFile = 1;
% show waveform metrics
p.Results.showFigures = 1;
% ISI threshold
p.ISIthreshold = 2/10^3; % unit:sec

% waveform import parameters
p.gwfparams.dataType = 'int16';            % Data type of .dat file (this should be BP filtered)
p.gwfparams.wfWin = [-40 41];              % Number of samples before and after spiketime to include in waveform
p.gwfparams.nWf = 2000;                    % Number of maximum waveforms per unit to pull out
p.gwfparams.nWf_plot = 1000;               % Number of maximum waveforms to plot
p.gwfparams.nCh = 385;


metaFile = dir(fullfile(apBinFolder,'*.ap.meta'));
ap_metaFileName = metaFile.name;
% Parse the corresponding metafile
apMetaFile = strcat(apBinFolder,'\',ap_metaFileName);
ap_meta = ReadMeta(apMetaFile);
Fs = SampRate(ap_meta);
% get digital to analog converter to uV
converter = ADconvert(meta);
p.gwfparams.ap_ADC = converter.AP(1);

% folderName is the folder where .ap.bin and ks4 folder is
% apTimeStamps should be the synchronized time with the NIDQ
% here I use raw time just as example, which is wrong in many case
nSamp = str2double(ap_meta.fileSizeBytes) / 2 / str2double(ap_meta.nSavedChans);
apTimeStamps = (1:length(nSamp))/Fs;
clusterInfo = clusterImport_ks4(folderName,Fs,apTimeStamps,p);

