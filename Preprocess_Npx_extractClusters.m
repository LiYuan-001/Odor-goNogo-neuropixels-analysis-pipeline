% This code is to extract good cluster timestamp, waveform and basic info
% after phy manual curation (only manual labled good in phy)
% Main goal is to assign spike to synchronized time
% After this step, anything related to good cluster should be reading ts
% file generated from this code
% Li Yuan, Oct-30-2022
function Preprocess_Npx_extractClusters(inFile,AnalyzeSes)

p.saveFile = 1;
p.savePlot = 1;
p.saveDir = 'N:\yuan\AnalysisResults'; % change it to your own folder
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
    
    if p.saveFile
        savedir = sprintf('%s%s%s%s%s',p.saveDir,'\ProcessedData\',sessInfo(sessInd).animalID,'-date-',sessInfo(sessInd).recDate); 
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
            savedir_plot = sprintf('%s%s%s%s%s%s%s',p.saveDir,'\Figures\',sessInfo(sessInd).animalID,'-date-',sessInfo(sessInd).recDate,'\SpikeShape\',probeName);
            if ~exist(savedir_plot, 'dir')
                mkdir(savedir_plot);
            end
            delete(strcat(savedir_plot,'\*'));
        end

        npx_Cluster.rat = sessInfo(sessInd).animalID;
        npx_Cluster.day = sessInfo(sessInd).recDate;
        npx_Cluster.DepthUnit = 'um';
        npx_Cluster.TsUnit = 'sec';

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
        p.savedir_plot = savedir_plot;
        clusterInfo = clusterImport_ks4_good_multi(ksPath,Fs,syncTime.Npx_timeStamps.(probeName).AP_Ts_Sync,p);
        clusterNum = length(clusterInfo);

        for k = 1:clusterNum
            clusterInfo{k}.avgRate = length(clusterInfo{k}.spkInd)./rec_Length;
        end
        
        npx_Cluster.clusterNum = clusterNum;
        npx_Cluster.clusterInfo = clusterInfo;
        
        save(sprintf('%s%s%s',savedir,'\npx_Cluster_',probeName),'npx_Cluster');
        close all
        clear npx_Cluster clusterInfo
        fprintf('Finished analysis for %s\n', probeName);
    end
    
    fprintf('Finished analysis for session %d\n',sessInd);
end