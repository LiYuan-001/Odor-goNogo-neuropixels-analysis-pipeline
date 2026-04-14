% This code is to extract good cluster timestamp, waveform and basic info
% after phy manual curation (only manual labled good in phy)
% Main goal is to assign spike to synchronized time
% After this step, anything related to good cluster should be reading ts
% file generated from this code
% Li Yuan, Oct-30-2022
function Preprocess_Npx_goodClusters(inFile,AnalyzeSes)

p.savePlot = 1;
p.writeToFile = 1;
% show waveform metrics
p.Results.showFigures = 1;
% ISI threshold
p.ISIthreshold = 2/10^3; % unit:sec

% waveform import parameters
p.gwfparams.dataType = 'int16';            % Data type of .dat file (this should be BP filtered)
p.gwfparams.wfWin = [-40 41];              % Number of samples before and after spiketime to include in waveform
p.gwfparams.nWf = 1000;                    % Number of maximum waveforms per unit to pull out
p.gwfparams.nWf_plot = 500;               % Number of maximum waveforms to plot
p.gwfparams.nCh = 385;
            
% Read in input information
sessInfo = SessInfoImport_Npx(inFile);

for sessInd = AnalyzeSes(1:end)
    close all
    
    if p.writeToFile
        savedir = strcat(sessInfo(sessInd).npx_path,'\Processed');
        if ~exist(savedir, 'dir')
            mkdir(savedir);
        end
    end

    % load corrected TIMESTAMPS
    syncFile = fullfile(savedir,'syncTime.mat');
    syncTime = load(syncFile);
    

    %% start processing each probe
    for m = 1:length(sessInfo(sessInd).probeID)
        probeName =  sprintf('imec%d',sessInfo(sessInd).probeID(m));
        if p.savePlot
            % directory for plot figures
            % generate a folder for each rat eah day under the current folder
            savedir_1 = sprintf('%s%s%s%s%s%s','N:\yuan\AnalysisResults\Figures\',sessInfo(sessInd).animalID,'-date-',sessInfo(sessInd).recDate,'\SpikeShape\',probeName);
            if ~exist(savedir_1, 'dir')
                mkdir(savedir_1);
            end
            delete(strcat(savedir_1,'\*'));
        end

        npx_goodCluster.rat = sessInfo(sessInd).animalID;
        npx_goodCluster.day = sessInfo(sessInd).recDate;
        npx_goodCluster.DepthUnit = 'um';
        npx_goodCluster.TsUnit = 'sec';

        % get digital to analog converter
        p.gwfparams.ap_ADC = syncTime.Npx_timeStamps.(probeName).converter.AP(1);
    
        % get the path for the catGT data and kilosort folder (ks4)
        % I always let the folder with name catGT*imecxx

        ksPathTemp = dir(fullfile(sessInfo(sessInd).npx_path,strcat('catgt*',probeName)));
        ksPath = strcat(ksPathTemp.folder,'\',ksPathTemp.name);
        metaFileTemp = dir(fullfile(ksPath,'*.ap.meta'));

        % Parse the corresponding metafile
        ap_meta = ReadMeta(strcat(metaFileTemp.folder,'\',metaFileTemp.name));
        Fs = SampRate(ap_meta);
        rec_Length = str2double(ap_meta.fileTimeSecs);        
        
        % combined multiple code together so it is a bit redundant
        % Li YUAN
        
        % read good clusters
        % use referenced corrected Ts
        p.sessInfo = sessInfo(sessInd);
        p.savedir_1 = savedir_1;
        clusterInfo = clusterImport_ks4(ksPath,Fs,syncTime.Npx_timeStamps.(probeName).AP_Ts_Sync,p);
        clusterNum = length(clusterInfo);

        for k = 1:clusterNum
            clusterInfo{k}.avgRate = length(clusterInfo{k}.spkInd)./rec_Length;
        end
        
        npx_goodCluster.clusterNum = clusterNum;
        npx_goodCluster.clusterInfo = clusterInfo;
        
        save(sprintf('%s%s%s',savedir,'\npx_goodCluster_',probeName),'npx_goodCluster');
        close all
        clear npx_goodCluster clusterInfo
        fprintf('Finished analysis for %s\n', probeName);
    end
    
    fprintf('Finished analysis for session %d\n',sessInd);
end