% function edited from kilosort
% return only good unit
% Li Yuan, UCSD, 2023-Jan-18
function wf = getWaveForms_3(gwfparams,p)
% function wf = getWaveForms(gwfparams)
%
% Extracts individual spike waveforms from the raw datafile, for multiple
% clusters. Returns the waveforms and their means within clusters.
%
% Contributed by C. Schoonover and A. Fink
%
% % EXAMPLE INPUT
% gwfparams.dataDir = '/path/to/data/';    % KiloSort/Phy output folder
% gwfparams.fileName = 'data.dat';         % .dat file containing the raw 
% gwfparams.dataType = 'int16';            % Data type of .dat file (this should be BP filtered)
% gwfparams.nCh = 32;                      % Number of channels that were streamed to disk in .dat file
% gwfparams.wfWin = [-40 41];              % Number of samples before and after spiketime to include in waveform
% gwfparams.nWf = 2000;                    % Number of waveforms per unit to pull out
% gwfparams.spikeTimes =    [2,3,5,7,8,9]; % Vector of cluster spike times (in samples) same length as .spikeClusters
% gwfparams.spikeClusters = [1,2,1,1,1,2]; % Vector of cluster IDs (Phy nomenclature)   same length as .spikeTimes
%
% % OUTPUT
% wf.unitIDs                               % [nClu,1]            List of cluster IDs; defines order used in all wf.* variables
% wf.spikeTimeKeeps                        % [nClu,nWf]          Which spike times were used for the waveforms
% wf.waveForms                             % [nClu,nWf,nCh,nSWf] Individual waveforms
% wf.waveFormsMean                         % [nClu,nCh,nSWf]     Average of all waveforms (per channel)
%                                          % nClu: number of different clusters in .spikeClusters
%                                          % nSWf: number of samples per waveform
%
% % USAGE
% wf = getWaveForms(gwfparams);

wfNSamples = length(p.gwfparams.wfWin(1):p.gwfparams.wfWin(end));

% spike time file
spk_TsIndFile = fullfile(gwfparams.dataDir,'spike_times.npy');
spk_tsInd = double(readNPY(spk_TsIndFile));
% load each spike's assigned cluster index
spk_ClusterIndFile = fullfile(gwfparams.dataDir,'spike_clusters.npy');
spk_clusterInd = readNPY(spk_ClusterIndFile);

% load cluster Info
cluInfo = tdfread(fullfile(gwfparams.dataDir,'cluster_info.tsv'));
unitIDs = cluInfo.cluster_id;
numUnits = size(unitIDs,1);
goodIDs = [];
ID_ch = [];
ID_depth = [];
for curUnitInd=1:numUnits
    curUnitID = unitIDs(curUnitInd);
    % only write down good units
    if strcmp(cluInfo.group(curUnitInd,:),'good ')
        goodIDs = [goodIDs;curUnitID];
        ID_ch = [ID_ch;cluInfo.ch(curUnitInd)];
        ID_depth = [ID_depth;cluInfo.depth(curUnitInd)];
    end
end

num_goodUnits = length(goodIDs);

% maxNum = min(p.gwfparams.nWf,length(spkInd));
spikeTimeKeeps = nan(1,p.gwfparams.nWf);
waveForms = nan(1,p.gwfparams.nWf,wfNSamples);
waveFormsMean = nan(1,wfNSamples);

for curUnitInd=1:num_goodUnits
    curUnitID = goodIDs(curUnitInd);
    curUnitCh = ID_ch(curUnitInd);
    % only write down good units
    curSpikeTimes = spk_tsInd(spk_clusterInd==curUnitID);
    curUnitnSpikes = size(curSpikeTimes,1);
    spikeTimesRP = curSpikeTimes(randperm(curUnitnSpikes));
    spikeTimeKeeps(curUnitInd,1:min([gwfparams.nWf curUnitnSpikes])) = sort(spikeTimesRP(1:min([gwfparams.nWf curUnitnSpikes])));
    for curSpikeTime = 1:min([gwfparams.nWf curUnitnSpikes])
        tmpWf = mmf.Data.x(1:gwfparams.nCh,spikeTimeKeeps(curUnitInd,curSpikeTime)+gwfparams.wfWin(1):spikeTimeKeeps(curUnitInd,curSpikeTime)+gwfparams.wfWin(end));
        waveForms(curUnitInd,curSpikeTime,:) = double(tmpWf(curUnitCh+1,:));
    end
    waveFormsMean(curUnitInd,:) = squeeze(nanmean(waveForms(curUnitInd,:,:),2));
%     disp(['Completed ' int2str(curUnitInd) ' units of ' int2str(num_goodUnits) '.']);
end


% Package in wf struct
wf.unitIDs = goodIDs;
wf.unitChan = ID_ch;
wf.unitDepth = ID_depth;
wf.spikeTimeKeeps = spikeTimeKeeps;
wf.waveForms = waveForms;
wf.waveFormsMean = waveFormsMean;
wf.spikeTimes = curSpikeTimes;

end