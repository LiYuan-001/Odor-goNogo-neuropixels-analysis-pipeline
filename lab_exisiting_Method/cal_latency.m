%% cal_latency

function [peak_posi,continous_latency,amp_latency]=cal_latency(psth_data,time_bin,cri_z,cri_bin,cri_amp,bin_size)
% psth_data should be z-scored
% time_bin, the range to get latency
% cri_bin, the number of bins exceeding the criteria
% cri_amp, amplitude to detect the latency, e.g., 30%
% bin_size pf PSTH

data=psth_data;
c_num=cri_bin;

% data=tmpS1;
% c_num=cri_bin2;
% time_bin=time_tmp;

%% peak latency
tmpd=data(time_bin);
[~,pM]=max(tmpd);

x=1:pM+4;
y=data(x+time_bin(1));
t_x=linspace(x(1),x(end),bin_size*length(x));% ms resolution
itp_y = interp1(x,y,t_x);
[~,peak_posi]=max(itp_y);



%% continous_latency
%data=[0 1 0 0 1 0 0 1 0 1 1 0 0 1 1 1 0 1 1 0 1 1 1 1 1];
%c_num2=10;

tmpd=data(time_bin)>cri_z;
for i=1:length(tmpd)-c_num+1
    x(i)=sum(tmpd(i:c_num-1+i));
end

continous_bin=find(x>=c_num,1,'first');
if ~isempty(continous_bin)
    x=1:continous_bin+4;
    y=data(x+time_bin(1));
    t_x=linspace(x(1),x(end),bin_size*length(x));% ms resolution
    itp_y = interp1(x,y,t_x);
    
    c_num2=10;% continous_10ms
    tmpd=itp_y>cri_z;
    for i=1:length(tmpd)-c_num2+1
        x(i)=sum(tmpd(i:c_num2-1+i));
    end
    
    continous_latency=find(x>cri_amp,1,'first');
    if isempty(continous_latency)
        continous_latency=NaN;
    end
else
    continous_latency=NaN;
end


%% amp_latency
if ~isnan(continous_latency)
    tmpd=data(time_bin);
    ndata=data/max(tmpd);
    tmpd=tmpd/max(tmpd);
    tmp_bin=find(tmpd>cri_amp,1,'first');
    
    x=1:tmp_bin+4;
    y=ndata(x+time_bin(1));
    t_x=linspace(x(1),x(end),bin_size*length(x));% ms resolution
    itp_y = interp1(x,y,t_x);
    tmp_bin=find(itp_y>cri_amp,1,'first');
    amp_latency=tmp_bin;
else
    amp_latency=NaN;
end

if isempty(amp_latency)
    amp_latency=NaN;
end

end

