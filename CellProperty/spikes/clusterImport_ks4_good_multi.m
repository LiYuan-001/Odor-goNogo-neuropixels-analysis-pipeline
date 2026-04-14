% load spike cluster index, time sample index, cluster label and waveform
% after manual phy curation
% keep information for good clusters only
% Li Yuan, Jul-22-2022, UCSD
% this version if for kilosrt result saved in child folder called
% 'kilosort4'
% read in manual labeled 'good' and 'multi' clusters
% Li Yuan, Oct-30-2025
function cluster = clusterImport_ks4_good_multi(filePath,Fs,TimeStamps,p)
   
if ~isfield(p,'savePlot')
    p.savePlot = 0;
end

% get file name for waveform raw file
gwfparams.dataDir = filePath;
oldFolder = cd;
cd(filePath);
gwfparams.fileName = dir('*.ap.bin');
% route back to original folder

chanMapName = dir('*ChanMap.mat');
chanMapTemp = load(chanMapName.name);

cd(oldFolder);
fileName = fullfile(gwfparams.dataDir,gwfparams.fileName.name);           
filenamestruct = dir(fileName);
dataTypeNBytes = numel(typecast(cast(0, p.gwfparams.dataType), 'uint8')); % determine number of bytes per sample
nSamp = filenamestruct.bytes/(p.gwfparams.nCh*dataTypeNBytes);  % Number of samples per channel


% specify file name
paramFile = fullfile(filePath, 'kilosort4\params.py');
spk_ClusterIndFile = fullfile(filePath,'kilosort4\spike_clusters.npy');
spk_TsIndFile = fullfile(filePath,'kilosort4\spike_times.npy');
if exist(fullfile(filePath, 'kilosort4\cluster_group.csv')) 
    cgsFile = fullfile(filePath, 'kilosort4\cluster_group.csv');
end
if exist(fullfile(filePath, 'kilosort4\cluster_group.tsv')) 
    cgsFile = fullfile(filePath, 'kilosort4\cluster_group.tsv');
end

chanMapFile = fullfile(filePath, 'kilosort4\channel_map.npy');
chanMap = readNPY(chanMapFile);

% load x,y coords for recording channels
coords = readNPY(fullfile(filePath, 'kilosort4\channel_positions.npy'));
xcoords = coords(:,1);
ycoords = coords(:,2); 
% load each spike's assigned cluster index
spk_clusterInd = readNPY(spk_ClusterIndFile);
% load each spike sample index
spk_tsInd = double(readNPY(spk_TsIndFile));
% load templates
temps = readNPY(fullfile(filePath, 'kilosort4\templates.npy'));
% load whitening mat inv
winv = readNPY(fullfile(filePath, 'kilosort4\whitening_mat_inv.npy'));
% load templates
spikeTemplates = readNPY(fullfile(filePath, 'kilosort4\spike_templates.npy')); % note: zero-indexed
% load scaling amps
tempScalingAmps = readNPY(fullfile(filePath, 'kilosort4\amplitudes.npy'));
[spikeAmps, ~, ~, ~, ~, ~, ~] = ...
    templatePositionsAmplitudes(temps, winv, ycoords, spikeTemplates, tempScalingAmps);


% load cluster label from ohy (good, MUA, noise)
% - 0 = noise
% - 1 = mua
% - 2 = good
% - 3 = unsorted
[cids, cgs] = readClusterGroupsCSV(cgsFile);
% load cluster Info
cluInfo = tdfread(fullfile(filePath,'kilosort4\cluster_info.tsv'));

% get good clusters spike ind, time etc
goodInd = find(cgs == 2);
% get multi clusters spike ind, time etc
multiInd = find(cgs == 1);
allInd = [goodInd(:);multiInd(:)];

for k = 1:length(allInd)
    indTemp = cids(allInd(k));
    cluster{k}.ID = indTemp;
    cluster{k}.manualLabel = cgs(allInd(k));
    % remove repeated ind potentially from merge
    index = find(spk_clusterInd == indTemp);
    if length(index) > length(unique(index))
        index = unique(index);
        sprintf('Repeated spikes removed')
    end
    cluster{k}.spkInd = spk_tsInd(index);
    cluster{k}.spkInd = cluster{k}.spkInd(cluster{k}.spkInd<=length(TimeStamps));
    cluster{k}.spkTs = TimeStamps(cluster{k}.spkInd);
    if size(cluster{k}.spkTs,2) > size(cluster{k}.spkTs,1)
        cluster{k}.spkTs = cluster{k}.spkTs';
    end
    cluster{k}.spkAmps = spikeAmps(index);
    % find the matching template for this cluster
    % if a cluster was merged from multiple templates, multiple templates
    % number show up
    % use the most spikes as an estimate
    templateIndTemp = spikeTemplates(index);
    % spikeTemplates start from 0
    templateIndTemp_2 = mode(templateIndTemp)+1;
%     cluster{i}.mass_depth = templateDepths(templateIndTemp_2);
%     cluster{i}.template_waveform = waveforms(templateIndTemp_2,:)*1000;
    % assign channels to cluster
    ind2 = (cluInfo.cluster_id == indTemp);
    % if ~strcmpi(cluInfo.group(ind2,1:4),'good')
    %     error('Not same cluster')
    % end
    cluster{k}.channel.channel_ID = cluInfo.ch(ind2);
    cluster{k}.channel.channel_depth = cluInfo.depth(ind2);
    cluster{k}.channel.xcoords = chanMapTemp.xcoords(cluInfo.ch(ind2)+1);
    cluster{k}.channel.ycoords = chanMapTemp.ycoords(cluInfo.ch(ind2)+1);
    
    if cluster{k}.channel.ycoords ~= cluster{k}.channel.channel_depth
        error('error in channel depth assignment')
    end
    
    % get raw waveforms   
    mmf = memmapfile(fileName, 'Format', {p.gwfparams.dataType, [p.gwfparams.nCh nSamp], 'x'});
    wf = getWaveForms3(mmf,cluster{k}.spkInd,cluInfo.ch(ind2),p);
    wf.waveForms = wf.waveForms * p.gwfparams.ap_ADC(1)*10^6;% unit uV
    wf.waveFormsMean = wf.waveFormsMean * p.gwfparams.ap_ADC(1)*10^6; % unit uV
    % calculate the waveform property etc    
    metrics = clusterMetrics(cluster{k}.spkTs,wf.waveFormsMean,Fs,p);
    cluster{k}.waveform = wf.waveFormsMean;
    cluster{k}.waveMetrics = metrics;
    
    h = figure(1);
    h.Position = [100,100,900,900];
    subplot(2,2,2);
    if size(wf.waveForms,1) > p.gwfparams.nWf_plot
        indTemp = randperm(size(wf.waveForms,1));
        waveIdx = indTemp(1:p.gwfparams.nWf_plot);
    else
        waveIdx = 1:size(wf.waveForms,1);
    end
    wavePlot = wf.waveForms(waveIdx,:);    
    plot(wavePlot');
    xlim([0 length(wf.waveFormsMean)])
    TITLE1 = sprintf('%s%d%s%d%s%d','Chan-',cluster{k}.channel.channel_ID,'-Depth-',cluster{k}.channel.ycoords,'-xCoord-',cluster{k}.channel.xcoords);
    TITLE2 = sprintf('%s%d%s','random ',length(waveIdx),' spikes');
    title({TITLE1;TITLE2});
    
    subplot(2,2,1)
    plot([p.gwfparams.wfWin(1):p.gwfparams.wfWin(end)]/Fs*1000,cluster{k}.waveform,'r','LineWidth',3);
    xlabel('ms')
    ylabel('uV')
    if isfield(p,'sessInfo')
        TITLE1 = sprintf('%s%s%s%s','Rat-',p.sessInfo.animalID,'-Day-',p.sessInfo.recDate);
    else
        TITLE1 = [];
    end
    TITLE2 = 'Average waveform';
    title({TITLE1;TITLE2});
    hold on
    waveMax = max(cluster{k}.waveform);
    waveRange = max(cluster{k}.waveform)-min(cluster{k}.waveform);
    if waveRange < 100
        textPos = 0:5:40;
    elseif waveRange < 200
        textPos = 0:15:90;
    else
        textPos = 0:30:180;
    end
    if cgs(allInd(k)) == 2
        text(0.5,waveMax-textPos(1),sprintf('%s','Good'))
    else
        text(0.5,waveMax-textPos(1),sprintf('%s','Multi'))
    end
    text(0.5,waveMax-textPos(2),sprintf('%s%1.3f','duration: ',cluster{k}.waveMetrics.duration*10^3))
    % text(0.5,waveMax-textPos(2),sprintf('%s%1.3f','halfWidth: ',cluster{k}.waveMetrics.halfWidth*10^3))
    text(0.5,waveMax-textPos(3),sprintf('%s%1.3f','pt_ratio: ',cluster{k}.waveMetrics.pt_ratio),'Interpreter','None')
    text(0.5,waveMax-textPos(4),sprintf('%s%3.1f','amp: ',cluster{k}.waveMetrics.amp))
    % text(0.5,waveMax-textPos(5),sprintf('%s%1.2f','contam\_Ratio: ',cluster{k}.waveMetrics.acg.contamination))
    text(0.5,waveMax-textPos(5),sprintf('%s%2.2f','isi\_violation: ',cluster{k}.waveMetrics.isi_violationPct))
    
                
    %plot wide ACG
    subplot(2,2,3)
    plot(cluster{k}.waveMetrics.acg.acg_wide,'k')
    title('ACG wide (1s, 1ms bins)')
    xlim([0 length(cluster{k}.waveMetrics.acg.acg_wide)])
    set(gca,'XTick',[0 500 1000],'XTickLabel',{'-500','0','500'});
    xlabel('ms')
    % plor acg narrow
    subplot(2,2,4)
    midSize = ceil(length(cluster{k}.waveMetrics.acg.acg_narrow)/2);
    plot(-100:0.5:99.5,cluster{k}.waveMetrics.acg.acg_narrow,'r')
    % xlim([0 length(cluster{k}.waveMetrics.acg.acg_narrow)])
    % set(gca,'XTick',[0 100 200],'XTickLabel',{'-100','0','100'});
    xlabel('ms')
    hold on
    
    a = cluster{k}.waveMetrics.acg_fit_params.acg_tau_decay;
    b = cluster{k}.waveMetrics.acg_fit_params.acg_tau_rise;
    c = cluster{k}.waveMetrics.acg_fit_params.acg_c;
    d = cluster{k}.waveMetrics.acg_fit_params.acg_d;
    e = cluster{k}.waveMetrics.acg_fit_params.acg_asymptote;
    f = cluster{k}.waveMetrics.acg_fit_params.acg_refrac;
    g = cluster{k}.waveMetrics.acg_fit_params.acg_tau_burst;
    h = cluster{k}.waveMetrics.acg_fit_params.acg_h;
    
    x_fit = 1:0.1:100;
    fiteqn = max(c*(exp(-(x_fit-f)/a)-d*exp(-(x_fit-f)/b))+h*exp(-(x_fit-f)/g)+e,0);
    line([-flip(x_fit),x_fit],[flip(fiteqn),fiteqn],'linewidth',1.5,'color',[0,0,0,0.7],'HitTest','off')
    
    title('ACG narrow (100ms, 0.5ms bins)')
    xlabel(sprintf('%s%1.2f%s%1.2f','ACG tau rise: ',cluster{k}.waveMetrics.acg_fit_params.acg_tau_rise, ' decay: ',cluster{k}.waveMetrics.acg_fit_params.acg_tau_decay));
    %                 xlim([0 length(clusterInfo{k}.acg_metrics.acg_narrow)])
    %                 set(gca,'XTick',[-50,0,50],'XTickLabel',{'-50','0','50'});
    
    if p.savePlot == 1
        figure(1)
        figName = sprintf('%s%s%s%s%s%s%d%s%d%s%d',p.savedir_plot,'\Rat-',p.sessInfo.animalID,'-Date-',p.sessInfo.recDate,'-ycoord-',cluster{k}.channel.channel_depth,'-xCoord-',cluster{k}.channel.xcoords,'-KsID-',cluster{k}.ID);
        print(figName,'-dpng','-r300');
    end
    close all    
    clear wf
    fprintf('%s%d\n','Finish cluster: ',k);
end

% load waveforms for all good clusters


end