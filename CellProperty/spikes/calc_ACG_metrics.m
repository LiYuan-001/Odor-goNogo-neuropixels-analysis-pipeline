function acg_metrics = calc_ACG_metrics(ts,Fs,p)
% Two autocorrelograms are calculated:  narrow (100ms, 0.5ms bins) and wide (1s, 1ms bins) using the CCG function (for speed)
%
% Further three metrics are derived from these:
%
% Theta modulation index:
%    Computed as the difference between the theta modulation trough (defined as mean of autocorrelogram bins 50-70 ms)
%    and the theta modulation peak (mean of autocorrelogram  bins 100-140ms) over their sum
%    Originally defined in Cacucci et al., JNeuro 2004
%
% BurstIndex_Doublets:
%    max bin count from 2.5-8ms normalized by the average number of spikes in the 8-11.5ms bins
%
% BurstIndex_Royer2012:
%    Burst index is determined by calculating the average number of spikes in the 3-5 ms bins of the spike
%    autocorrelogram divided by the average number of spikes in the 200-300 ms bins.
%    Metrics introduced in Royer et al. Nature Neuroscience 2012, and adjusted in Senzai & Buzsaki, Neuron 2017.

% By Peter Petersen
% petersen.peter@gmail.com
% Last edited: 06-10-2020

bins_wide = 500;
acg_wide = zeros(bins_wide*2+1,1);
bins_narrow = 100;
acg_narrow = zeros(bins_narrow*2+1,1);
% disp('Calculating narrow ACGs (100ms, 0.5ms bins) and wide ACGs (1s, 1ms bins)')
% tic
% ccg REQUIRES MEX FILE 
acg_wide = ACG_2(ts,'binSize',0.001,'duration',0.5,'norm','rate');
acg_narrow = ACG_2(ts,'binSize',0.0005,'duration',0.1,'norm','rate');
% Metrics from narrow
BurstIndex_Doublets = max(acg_narrow(bins_narrow+1+5:bins_narrow+1+16))/mean(acg_narrow(bins_narrow+1+16:bins_narrow+1+23));
% Metrics from wide
ThetaModulationIndex = (mean(acg_wide(bins_wide+1+100:bins_wide+1+140)) - mean(acg_wide(bins_wide+1+50:bins_wide+1+70)))/(mean(acg_wide(bins_wide+1+50:bins_wide+1+70))+mean(acg_wide(bins_wide+1+100:bins_wide+1+140)));
BurstIndex_Royer2012 = mean(acg_wide(bins_wide+1+3:bins_wide+1+5))/mean(acg_wide(bins_wide+1+200:bins_wide+1+300));
% contamination rate
violationBin = round(10^3*p.ISIthreshold/1);
contamination = mean(acg_wide(bins_wide+1:bins_wide+1+violationBin))/mean(acg_wide(bins_wide+1:end));
% toc

% % log 10 ACG
% 
% intervals = -3:0.04:1;
% intervals2 = intervals(1:end-1)+.02;
% acg.log10 = zeros(length(intervals2),1);
% acg.log10_bins = 10.^intervals2';
% acg_log10 = zeros(length(intervals2),1);
% 
% parallel_toolbox_installed = isToolboxInstalled('Parallel Computing Toolbox'); % Validating that Parallel Computing Toolbox is installed
% parallel_toolbox_installed = 0;
% if parallel_toolbox_installed
%     gcp;
% 
%         ACGlog = zeros(1,length(intervals)-1);
%         i = 1;
%         test = 1;
%         while test > 0
%             ISIs = log10(ts(i+1:end)-ts(1:end-i));
%             [N,~] = histcounts(ISIs,intervals);
%             ACGlog = ACGlog+N;
%             i = i+1;
%             test = any(ISIs<intervals(end));
%         end
%         acg_log10 = ACGlog./(diff(10.^intervals))/length(ts);
% else
% 
%         ACGlog = zeros(1,length(intervals)-1);
%         i = 1;
%         test = 1;
%         while test > 0
%             ISIs = log10(ts(i+1:end)-ts(1:end-i));
%             [N,~] = histcounts(ISIs,intervals);
%             ACGlog = ACGlog+N;
%             i = i+1;
%             test = any(ISIs<intervals(end));
%         end
%         acg_log10 = ACGlog./(diff(10.^intervals))/length(ts);
% 
% end    
    
acg_metrics.acg_wide = acg_wide;
acg_metrics.acg_narrow = acg_narrow;
% acg_metrics.acg_log10 = acg_log10;
acg_metrics.contamination = contamination;
acg_metrics.thetaModulationIndex = ThetaModulationIndex;
acg_metrics.burstIndex_Royer2012 = BurstIndex_Royer2012;
acg_metrics.burstIndex_Doublets = BurstIndex_Doublets;

