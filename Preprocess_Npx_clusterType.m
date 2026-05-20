% This code is to read in extracted clusters and decide whether this is a
% good cluster for the rest analysis 
% I separate it from 'Preprocess_Npx_extractClusters.m' because I might set
% different threshold, but I do not want to go through extract spikes each
% time
% Li YUAN, Tohoku, 2026-Apr-20
function Preprocess_Npx_clusterType(inFile,AnalyzeSes)

p.saveFile = 1;
p.savePlot = 1;
p.saveDir = 'N:\yuan\AnalysisResults'; % change it to your own folder

% ISI threshold
p.isi_Thres = 2/10^3; % unit:sec
p.isi_vio_Thres = 2; % percent
p.spkWidth_Thres = 0.4*10^-3; % second
p.isiWindow = 0.05; % second
p.avgRate_Thres = 0.01; % Hz
          
% Read in input information
sessInfo = SessInfoImport_Npx(inFile);

for sessInd = AnalyzeSes(1:end)
    close all
    
    if p.saveFile
        savedir = sprintf('%s%s%s%s%s',p.saveDir,'\ProcessedData\',sessInfo(sessInd).animalID,'-date-',sessInfo(sessInd).recDate); 
        % if ~exist(savedir, 'dir')
        %     mkdir(savedir);
        % end
    end

    npx_Cluster_type.rat = sessInfo(sessInd).animalID;
    npx_Cluster_type.day = sessInfo(sessInd).recDate;
    npx_Cluster_type.DepthUnit = 'um';
    npx_Cluster_type.p = p;
    %% start processing each probe
    for m = 1:length(sessInfo(sessInd).probeID)
        probeName =  sprintf('imec%d',sessInfo(sessInd).probeID(m));

        if p.savePlot
            % directory for plot figures
            % generate a folder for each rat eah day under the current folder
            savedir_plot = sprintf('%s%s%s%s%s%s%s',p.saveDir,'\Figures\',sessInfo(sessInd).animalID,'-date-',sessInfo(sessInd).recDate,'\ClusterType\',probeName);
            if ~exist(savedir_plot, 'dir')
                mkdir(savedir_plot);
            end
            delete(strcat(savedir_plot,'\*'));
        end
           

        % load cluster spike time and waveform
        spikeFile = fullfile(savedir,sprintf('%s%s%s','npx_Cluster_',probeName,'.mat'));
        load(spikeFile);

        samp_ts = (npx_Cluster.gwfparams.wfWin(1):npx_Cluster.gwfparams.wfWin(2))./30000;
        clusterNum = npx_Cluster.clusterNum;
        xcoord = nan(clusterNum,1);
        ycoord = nan(clusterNum,1);
        cluster_ksID = nan(clusterNum,1);
        width = nan(clusterNum,1);
        isi_vio = nan(clusterNum,1);
        avgRate = nan(clusterNum,1);

        fprintf('%s %d clusters\n',probeName, clusterNum);
        for k = 1:clusterNum
            cluster_ksID(k) = npx_Cluster.clusterInfo{k}.ID;
            waveform = npx_Cluster.clusterInfo{k}.waveform;
            % tau_riseTemp = npx_Cluster.clusterInfo{k}.waveMetrics.acg_fit_params.acg_tau_rise;
            % tau_rise(k) = tau_riseTemp;
            width(k) = npx_Cluster.clusterInfo{k}.waveMetrics.duration;
            xcoord(k) = npx_Cluster.clusterInfo{k}.channel.xcoords;
            ycoord(k) = npx_Cluster.clusterInfo{k}.channel.ycoords;
            avgRate(k) = npx_Cluster.clusterInfo{k}.avgRate;

            tsp = npx_Cluster.clusterInfo{k}.spkTs;
            [violation_Pct,~] = isi_violations(tsp,p.isi_Thres, 0);
            isi_vio(k) = violation_Pct;

            h = figure;
            h.Position = [100,100,900,500];

            subplot(1,2,1)
            plot(samp_ts*10^3,waveform,'k','LineWidth',2);
            xlabel('ms')
            ylabel('uV')
            TITLE1 = sprintf('%s%d%s%d%s%d%s%d','Chan-',npx_Cluster.clusterInfo{k}.channel.channel_ID,'-ycoord-',ycoord(k),'-xCoord-',xcoord(k),'-ksID-',cluster_ksID(k));
            if width(k) >= p.spkWidth_Thres
                TITLE2 = sprintf('%s%1.3f%s','Wide. Duration: ',width(k)*10^3,' ms');
            else
                TITLE2 = sprintf('%s%1.3f%s','Narrow. Duration: ',width(k)*10^3,' ms');
            end
            title({TITLE1;TITLE2},'Interpreter','None');

            subplot(1,2,2)
            isi = diff(tsp)*10^3;   % ms

            binEdges = (0:p.isi_Thres:p.isiWindow)*10^3;
            histogram(isi, binEdges, ...
                'FaceColor', 'r', ...
                'EdgeColor', 'r');

            xline(p.isi_Thres*10^3,'k--')
            xlim([0 p.isiWindow]*10^3)
            xlabel('Inter-spike interval (ms)')
            ylabel('Count')

            TITLE1 = sprintf('%s%2.2f%s', 'ISI vio: ', isi_vio(k), '%');
            TITLE2 = sprintf('%s%3.2f%s', 'AvgRate: ', avgRate(k), ' Hz');
            title({TITLE1; TITLE2}, 'Interpreter', 'none');

            if p.savePlot == 1
                figName = sprintf('%s%s%s%s%s%s%d%s%d%s%d',savedir_plot,'\clusterType-',sessInfo(sessInd).animalID,'-Date-',sessInfo(sessInd).recDate,'-ycoord-',ycoord(k),'-xCoord-',xcoord(k),'-ksID-',cluster_ksID(k));
                print(figName,'-dpng','-r200');
            end
            close(h);
        end

        npx_Cluster_type.(probeName).xcoord = xcoord;
        npx_Cluster_type.(probeName).ycoord = ycoord;
        npx_Cluster_type.(probeName).cluster_ksID = cluster_ksID;
        npx_Cluster_type.(probeName).avgRate = avgRate;
        npx_Cluster_type.(probeName).width = width;
        npx_Cluster_type.(probeName).isi_vio = isi_vio;
        
        npx_Cluster_type.(probeName).wide_label = width >= p.spkWidth_Thres;
        npx_Cluster_type.(probeName).single_label = isi_vio <= p.isi_vio_Thres;

        clear npx_Cluster
        fprintf('Finished analysis for %s\n', probeName);
    end

    save(sprintf('%s%s%s',savedir,'\npx_Cluster_type'),'npx_Cluster_type');
    clear npx_Cluster_type
    close all
    fprintf('Finished analysis for session %d\n',sessInd);
end