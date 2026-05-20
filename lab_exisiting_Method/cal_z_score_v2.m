%% cal_z_score v2
%Modified 230414 JL
%%Changes condition from if std(data(base))==0 to <0.01 to avoid bug where
%%near-zero baseline yields huge z-scores (z=100 to 1000!)

function [z_data]= cal_z_score_v2(psth,bin_baseline)

base=bin_baseline;
data=psth;

if sum(data)==0
    z_data=zeros(1,length(psth));
else
    if std(data(base)) < 0.05 %This condition modified 230508
        z_data=(data-mean(data(base)))/std(data);
    elseif std(data)==0
        z_data=zeros(1,length(psth));
    else
        z_data=(data-mean(data(base)))/std(data(base));
    end
end

end




