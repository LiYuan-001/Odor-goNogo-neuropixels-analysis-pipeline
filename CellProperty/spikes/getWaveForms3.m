% function edited from kilosort
% return only good unit
% Li Yuan, UCSD, 2023-Jan-18
function wf = getWaveForms3(mmf,spkInd,ch,p)
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
maxNum = min(p.gwfparams.nWf,length(spkInd));
% spikeTimeKeeps = nan(1,maxNum);
waveForms = nan(maxNum,wfNSamples);
% waveFormsMean = nan(1,wfNSamples);

% decide to take out random spike index out
curUnitnSpikes = size(spkInd,1);
spikeTimesRP = spkInd(randperm(curUnitnSpikes));
spikeTimeKeeps = sort(spikeTimesRP(1:maxNum));
    
for curSpikeTime = 1:maxNum
    tmpWf = mmf.Data.x(1:p.gwfparams.nCh,spikeTimeKeeps(curSpikeTime)+p.gwfparams.wfWin(1):spikeTimeKeeps(curSpikeTime)+p.gwfparams.wfWin(end));
    waveForms(curSpikeTime,:) = double(tmpWf(ch+1,:));
end
waveFormsMean = nanmean(waveForms,1);

wf.waveForms = waveForms;
wf.waveFormsMean = waveFormsMean;

end