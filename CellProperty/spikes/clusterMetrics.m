function metrics = clusterMetrics(spkTs,waveFormsMean,Fs,p)
waveMeasure = waveformMetrics(waveFormsMean,Fs);
metrics.duration = waveMeasure.duration;
metrics.halfWidth = waveMeasure.halfWidth;
metrics.pt_ratio = waveMeasure.pt_ratio;
metrics.amp = waveMeasure.amp;

% calculate isi violation
[metrics.isi_violationPct,metrics.isi_violationRate] = isi_violations(spkTs,p.ISIthreshold, 0);

% calculate ACG
% save time
spkTs2 = spkTs(1:(min(length(spkTs),10^7)));
acg_metrics = calc_ACG_metrics(spkTs2,Fs,p);
metrics.acg = acg_metrics;
% calculate ACG fit
fit_params_out = fit_ACG(acg_metrics.acg_narrow);
metrics.acg_fit_params = fit_params_out;

end

function waveMeasure = waveformMetrics(wave,Fs)
peakTemp = max(wave);
troughTemp = min(wave);
if abs(peakTemp) <= abs(troughTemp)
% detect trough first and then peak
    [trough,trough_idx] = min(wave);
    [peak,peak_idx] = max(wave(trough_idx:end));
    peak_idx = peak_idx + trough_idx-1;
    waveMeasure.duration = (peak_idx - trough_idx)/Fs;
    
    threshold = wave(peak_idx)*0.5;
    thresh_idx1 = wave(peak_idx:end)<= threshold;
    thresh_idx1 = thresh_idx1(1) + peak_idx;
    thresh_idx2 = wave(1:peak_idx)<= threshold;
    thresh_idx2 = thresh_idx2(end);
    waveMeasure.halfWidth = (thresh_idx1 - thresh_idx2)/Fs;
    
else
    [peak,peak_idx] = max(wave);
    [trough,trough_idx] = min(wave(peak_idx:end));
    trough_idx = trough_idx + peak_idx-1;
    waveMeasure.duration = (trough_idx - peak_idx)/Fs;
    
    threshold = wave(trough_idx)*0.5;
    thresh_idx1 = wave(trough_idx:end)>= threshold;
    thresh_idx1 = thresh_idx1(1) + trough_idx;
    thresh_idx2 = wave(1:trough_idx)>= threshold;
    thresh_idx2 = thresh_idx2(end);
    waveMeasure.halfWidth = (thresh_idx1 - thresh_idx2)/Fs;
end

waveMeasure.pt_ratio = abs(wave(peak_idx)/wave(trough_idx));
waveMeasure.amp = abs(wave(peak_idx)-wave(trough_idx));


end


function [violation_Pct,violation_Rate] = isi_violations(spkTs,isi_threshold, min_isi)

% spike_train : array of spike times
%  min_time : minimum time for potential spikes
% max_time : maximum time for potential spikes
% isi_threshold : threshold for isi violation
% min_isi : threshold for duplicate spikes

if size(spkTs,2) > size(spkTs,1)
    spkTs = spkTs';
end
duplicate_spikes = diff(spkTs) <= min_isi;
duplicate_spikes2 = [0;duplicate_spikes];
spkTs_2 = spkTs(~duplicate_spikes2);

isi = diff(spkTs_2);
num_spikes = length(spkTs);
num_violations = sum(isi < isi_threshold);
violation_time = 2 * num_spikes * (isi_threshold - min_isi);

total_rate = length(spkTs)./(max(spkTs)-min(spkTs));
violation_Pct = 100*num_violations/num_spikes;
violation_rate = num_violations / violation_time;
violation_Rate = violation_rate / total_rate;

end
