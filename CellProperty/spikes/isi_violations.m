function [violation_Pct,violation_Rate] = isi_violations(spkTs,isi_threshold, min_isi)

% spike_train : array of spike times
% min_time : minimum time for potential spikes
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