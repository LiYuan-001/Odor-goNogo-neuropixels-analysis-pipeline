%% cal_peak

function [peak,peak_pm_bins]=cal_peak(psth_data,time_bin,bin_size,bin_range)
% psth_data should be spikes/sec data
% time_bin, the range to get peak
% bin_size pf PSTH

data=psth_data;
bins=round(bin_range/bin_size/2);

%data=tmpS2;
% time_bin=time_tmp;

%% peak
tmpd=data(time_bin);
[peak,pM]=max(tmpd);

%%  peak +/- bin
if pM-bins<0
    x=1:bins*2;
    peak_pm_bins=mean(tmpd(x));
else
    x=pM+time_bin(1)-bins:pM+time_bin(1)+bins-1;
    peak_pm_bins=mean(data(x));
end


end

