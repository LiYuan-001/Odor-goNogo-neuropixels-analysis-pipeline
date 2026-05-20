%%  gng_cal_history_V4

%230414 JL
%Based on gng_cal_history_V3 (function nested within gngSessionSum3)

%Part 1:
%Updates individual spike files previously calculated with cal_z_score
%Uses new cal_z_score_v2 to update history.ABCD_PSTH_gs_z and history.AB_PSTH_gs_z
%Also updates variables calculated from PSTH_gs_z e.g. history.his_PSTH_gs_z, history.AB_lat_peak_odor

%Part 2: Creates new history struct fields for response periods with new time windows
%E.g. Currently we have history.ABCD_signrankOdor or ABCD_signrankDelay which are calculated from 0-1s or 1-3s.
%We want to check some variations such as 0.5-1.5s, 0.5-2s, 0.5-3s
%To be named e.g. history.ABCD_signrank05to15, _05to2, _05to3



 
 
%%% config
bin=0.05;%sec, bin for psth
bin_size=50;%ms, bin for psth
win=[-2 6];%sec, bin for raster and psth(8 sec, 161 bins data)
spike=spike_time;%sec, time stamp
trgrange=10;% trial numbers for comparison between pre and post novel odor presentation

bin_baseline=21:40;% to cal z-score, -1 to 0 s
bin_odor=41:60;% 0 to 1 s
bin_odor2=51:70;% 0.5 to 1.5 s
bin_delay=81:100;%2 - 3
bin_delay2=71:100;%1.5 - 3
bin_rwd=101:120;%3 to 4s

%New bin ranges 230414
bin_odor05_2 = 51:80;
bin_odor05_3 = 51:100;

%cri_z=1.96;% p<0.05
cri_z=2.58;% p<0.01
%cri_z=3.25;% p<0.001

cri_bin1=3;% 150ms
cri_bin2=6;% 300ms

cri_amp=0.3;%30% for latency detection
bin_range=300;% peak +/- bin_range Soma et al., 2019


%% [1] History analysis - comparison between pre and post novel odor presentation
%%% cal raster and psth
history.his_name{1,1}='A pre';
history.his_name{1,2}='A post';
history.his_name{2,1}='B pre';
history.his_name{2,2}='B post';
for s=1:2% odor A and B
    oid=s;
    
    %%% get target time stamp
    trgid=history.time_stamp(:,3)==oid;
    tmptrg=history.time_stamp(trgid,:);%sec
    tmptrg=tmptrg(~isnan(tmptrg(:,7)),:);% removing immature trials
    tmptrg=tmptrg(tmptrg(:,8)==0,:);% removing manual rwd trials
    
    for ss=1:2% pre- and post- novel odor presetnation
        if ss==1% pre
            trgid=find(tmptrg(:,1)<0,trgrange,'last');
        elseif ss==2% post
            trgid=find(tmptrg(:,1)>0,trgrange,'first');
        end
        trg=tmptrg(trgid,2);%sec
        
        % --------------------------------------------------
        %%% cal rstar
        rst=ana_raster_ss(trg,win,spike);
        
        %%% cal psth
        psth = ana_peth_ss(rst,win,bin);
        
        %%% cal somooth psth
        bin_gs=50;%ms
        psth_gs=smoothing_gaussian_ss(psth,100,bin_gs,2);
        
        %%% cal z scored psth
        psth_gs_z=cal_z_score_v2(psth_gs,bin_baseline);
        % --------------------------------------------------
        history.his_Rater{s,ss}=rst;
        history.his_PSTH{s,ss}=psth(:);
        history.his_PSTH_gs{s,ss}=psth_gs(:);
        history.his_PSTH_gs_z{s,ss}=psth_gs_z(:);

        % --------------------------------------------------
        %%% detection of significant responses & latency detection
        tmpS1=history.his_PSTH_gs_z{s,ss};
        
        time_tmp=bin_odor;
        tmpd=tmpS1(time_tmp)>cri_z;
        [~,sig_id]=cal_continuous_num(tmpd,cri_bin1);
        tmp_sig1(s,ss)=sig_id;
        
        
        time_tmp=bin_delay;
        tmpd=tmpS1(time_tmp)>cri_z;
        [~,sig_id]=cal_continuous_num(tmpd,cri_bin2);
        tmp_sig2(s,ss)=sig_id;
        
        time_tmp=bin_rwd;
        tmpd=tmpS1(time_tmp)>cri_z;
        [~,sig_id]=cal_continuous_num(tmpd,cri_bin1);
        tmp_sig3(s,ss)=sig_id;
        
        %New bin ranges added 230414
        time_tmp=bin_odor2; %0.5-1.5s
        tmpd=tmpS1(time_tmp)>cri_z;
        [~,sig_id]=cal_continuous_num(tmpd,cri_bin1);
        tmp_sig4(s,ss)=sig_id;
        
        time_tmp=bin_odor05_2; %0.5-2s
        tmpd=tmpS1(time_tmp)>cri_z;
        [~,sig_id]=cal_continuous_num(tmpd,cri_bin1);
        tmp_sig5(s,ss)=sig_id;
        
        time_tmp=bin_delay2; %1.5-3s
        tmpd=tmpS1(time_tmp)>cri_z;
        [~,sig_id]=cal_continuous_num(tmpd,cri_bin1);
        tmp_sig6(s,ss)=sig_id;
        
        time_tmp=bin_odor05_3; %0.5-3s
        tmpd=tmpS1(time_tmp)>cri_z;
        [~,sig_id]=cal_continuous_num(tmpd,cri_bin1);
        tmp_sig7(s,ss)=sig_id;
        
        % --------------------------------------------------
    end
end
history.his_SigID_odor_Pre  =tmp_sig1(:,1);
history.his_SigID_delay_Pre =tmp_sig2(:,1);
history.his_SigID_rwd_Pre   =tmp_sig3(:,1);
history.his_SigID_odor_Post =tmp_sig1(:,2);
history.his_SigID_delay_Post=tmp_sig2(:,2);
history.his_SigID_rwd_Post  =tmp_sig3(:,2);

%Added 230414
history.his_SigID_odor05_15_Pre  =tmp_sig4(:,1);
history.his_SigID_odor05_2_Pre =tmp_sig5(:,1);
history.his_SigID_delay15_3_Pre   =tmp_sig6(:,1);
history.his_SigID_odor05_3_Pre =tmp_sig7(:,1);
history.his_SigID_odor05_15_Post =tmp_sig4(:,2);
history.his_SigID_odor05_2_Post=tmp_sig5(:,2);
history.his_SigID_delay15_3_Post  =tmp_sig6(:,2);
history.his_SigID_odor05_3_Post  =tmp_sig7(:,2);

%%230414- COMMENTED OUT history boundaries analysis below (redundant)
%%231205- BROUGHT BACK history boundaries analysis below - why did we
%%previously comment it out???
%% cal_history_data2 - Transferred from cal_history_data2.m 200711 Kei
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% cal_history_data2

%%% config
wind=500;%ms
spike=spike_time;
wins=[-1 -0.5;-0.5 0;0 0.5; 0.5 1;1 1.5; 1.5 2;2 2.5;2.5 3;3 3.5;3.5 4;4 4.5;4.5 5;5 5.5;5.5 6];
history.time_stamp_name{9,1}='Baseline Pre(Hz)';
history.time_stamp_name{10,1}='Baseline Pre(Hz)';
history.time_stamp_name{11,1}='Odor 0-0.5(Hz)';
history.time_stamp_name{12,1}='Odor 0.5-1(Hz)';
history.time_stamp_name{13,1}='Delay 1-1.5(Hz)';
history.time_stamp_name{14,1}='Delay 1.5-2(Hz)';
history.time_stamp_name{15,1}='Delay 2-2.5(Hz)';
history.time_stamp_name{16,1}='Delay 2.5-3(Hz)';
history.time_stamp_name{17,1}='Resp 3-3.5(Hz)';
history.time_stamp_name{18,1}='Resp 3.5-4(Hz)';
history.time_stamp_name{19,1}='Resp 4-4.5(Hz)';
history.time_stamp_name{20,1}='Resp 4.5-5(Hz)';
history.time_stamp_name{21,1}='Baseline Post(Hz)';
history.time_stamp_name{22,1}='Baseline Post(Hz)';


num_add_param=14;
history.time_stamp=history.time_stamp(:,1:8);

%% [1] History analysis - cal # of spike

tmpNumSpikes=zeros(size(history.time_stamp,1),size(wins,1));
for s=1:size(history.time_stamp,1)% # of trials
    
    trg=history.time_stamp(s,2);%sec get time stamp
    
    for ss=1:size(wins,1)% get number of spikes in the specific time period
        
        % --------------------------------------------------
        %%% cal rstar
        rst=ana_raster_ss(trg,wins(ss,:),spike);
        tmpNumSpikes(s,ss)=length(rst(~isnan(rst(:,1)) & rst(:,1)~=0,1))/wind*1000;%Hz
        
        % --------------------------------------------------
        
    end
end

history.time_stamp_DM=[history.time_stamp tmpNumSpikes];



%% [2] History analysis - categorizing trials pattern
%%% config
tmp_data_matrix=history.time_stamp_DM;
trm_range=10;
trm_range2=100;

%% [2-1] Boundary: First novel odor presentation
Boundary1=tmp_data_matrix(:,1)<0; %Pre-Boundary
Boundary2=tmp_data_matrix(:,1)>=0;%Post-Boundary


if exist('numOdor','var')
    
    if numOdor.ABCDsession ==4
        sRange = 1:4;
    elseif numOdor.ABCDsession ==2
        sRange = 3:4;
    end
else
    
    sRange = 1:4;
end


for s=sRange% ABCD
% for s=1:4%A to D


    tmpDM1=tmp_data_matrix(Boundary1 & tmp_data_matrix(:,3)==s,9:9+num_add_param-1);
    tmpDM2=tmp_data_matrix(Boundary2 & tmp_data_matrix(:,3)==s,9:9+num_add_param-1);
    
    if ~isempty(tmpDM1) && size(tmpDM1,1)-trm_range+1>0
        tmpDMpre=tmpDM1(size(tmpDM1,1)-trm_range+1:size(tmpDM1,1),:);
    elseif ~isempty(tmpDM1)
        tmpDMpre=[nan(trm_range-size(tmpDM1,1),num_add_param);tmpDM1];
    else
        tmpDMpre=nan(trm_range,num_add_param);
    end
    
    if size(tmpDM2,1)-trm_range2<0
        tmpDMpost=[tmpDM2(1:end,:);nan(trm_range2-size(tmpDM2,1),num_add_param)];
    else
        tmpDMpost=tmpDM2(1:trm_range2,:);
    end
    
    history.BDRY_First_novel{s}=[tmpDMpre;tmpDMpost];
end


%% [2-2] Boundary: First Rwd
Target_trial=find(history.time_stamp(:,3)==3 &  history.time_stamp(:,7)==1,1,'first');
Boundary1=[ones(Target_trial-1,1);zeros(size(tmp_data_matrix,1)-Target_trial+1,1)];%Pre-Boundary
Boundary2=[zeros(Target_trial-1,1);ones(size(tmp_data_matrix,1)-Target_trial+1,1)];%Post-Boundary

if isempty(Boundary1)
    for s=sRange %A to D
        history.BDRY_First_Rwd{s}=nan(trm_range+trm_range2,num_add_param);
    end
else
    for s=sRange %A to D
        tmpDM1=tmp_data_matrix(Boundary1 & history.time_stamp(:,3)==s,9:9+num_add_param-1);
        tmpDM2=tmp_data_matrix(Boundary2 & history.time_stamp(:,3)==s,9:9+num_add_param-1);
        
        if ~isempty(tmpDM1) && size(tmpDM1,1)-trm_range+1>0
            tmpDMpre=tmpDM1(size(tmpDM1,1)-trm_range+1:size(tmpDM1,1),:);
        elseif ~isempty(tmpDM1)
            tmpDMpre=[nan(trm_range-size(tmpDM1,1),num_add_param);tmpDM1];
        else
            tmpDMpre=nan(trm_range,num_add_param);
        end
        
        if ~isempty(tmpDM2) && size(tmpDM2,1)-trm_range2+1>0
            tmpDMpost=tmpDM2(1:trm_range2,:);
        elseif ~isempty(tmpDM2)
            tmpDMpost=[nan(trm_range2-size(tmpDM2,1),num_add_param);tmpDM2];
        else
            tmpDMpost=nan(trm_range2,size(tmpDM2,2));
        end
        
        history.BDRY_First_Rwd{s}=[tmpDMpre;tmpDMpost];
    end
end


%% [2-2] Boundary: First Qui
Target_trial=find(history.time_stamp(:,3)==4 &  history.time_stamp(:,7)==0,1,'first');
Boundary1=[ones(Target_trial-1,1);zeros(size(tmp_data_matrix,1)-Target_trial+1,1)];%Pre-Boundary
Boundary2=[zeros(Target_trial-1,1);ones(size(tmp_data_matrix,1)-Target_trial+1,1)];%Post-Boundary

if isempty(Boundary1)
    for s=1:4%A to D
        history.BDRY_First_Qui{s}=nan(trm_range+trm_range2,num_add_param);
    end
else
    for s=sRange %A to D
        tmpDM1=tmp_data_matrix(Boundary1 & history.time_stamp(:,3)==s,9:9+num_add_param-1);
        tmpDM2=tmp_data_matrix(Boundary2 & history.time_stamp(:,3)==s,9:9+num_add_param-1);
        
        if ~isempty(tmpDM1) && size(tmpDM1,1)-trm_range+1>0
            tmpDMpre=tmpDM1(size(tmpDM1,1)-trm_range+1:size(tmpDM1,1),:);
        elseif ~isempty(tmpDM1)
            tmpDMpre=[nan(trm_range-size(tmpDM1,1),num_add_param);tmpDM1];
        else
            tmpDMpre=nan(trm_range,num_add_param);
        end
        
        if ~isempty(tmpDM2) && size(tmpDM2,1)-trm_range2+1>0
            tmpDMpost=tmpDM2(1:trm_range2,:);
        elseif ~isempty(tmpDM2)
            tmpDMpost=[nan(trm_range2-size(tmpDM2,1),num_add_param);tmpDM2];
        else
            tmpDMpost=nan(trm_range2,size(tmpDM2,2));
        end
        
        history.BDRY_First_Qui{s}=[tmpDMpre;tmpDMpost];
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% END cal_history_data2.m

%% [2] General analysis - PSTH before novel odor presentation (AB block)
%%% config
trgrange=10;
% trial numbers for comparison between first and last X trials

%%% cal raster and psth
history.AB_name{1,1}='A all';
history.AB_name{1,2}=['A First',num2str(trgrange)];
history.AB_name{1,3}=['A Lastt',num2str(trgrange)];
history.AB_name{2,1}='B all';
history.AB_name{2,2}=['B First',num2str(trgrange)];
history.AB_name{2,3}=['B Lastt',num2str(trgrange)];

clear tmp_sig1 tmp_sig2 tmp_sig3 tmp_sig4 tmp_sig5 tmp_sig6 tmp_sig7
for s=1:2% odor A and B
    oid=s;
    
    %%% get target time stamp
    trgid=history.time_stamp(:,3)==oid; %Finds trials of either odor type A or B
    tmptrg=history.time_stamp(trgid,:);%sec
    tmptrg=tmptrg(~isnan(tmptrg(:,7)),:);% removing immature trials
    tmptrg=tmptrg(tmptrg(:,8)==0,:);% removing manual rwd trials
    
    AB_border=find(diff(tmptrg(:,4))<0); %<0 means trials that came before first novel odor; see history.time_stamp_name
    if isempty(AB_border)
        AB_border = size(tmptrg, 1);
    end
    tmptrg=tmptrg(1:AB_border,:);
    
    for ss=1:3% all, First, and Last
        if ss==1% all
            trgid=1:AB_border;
        elseif ss==2% first
            if AB_border>tmptrg
            trgid=1:trgrange;
            else
                trgid=1:round(AB_border/2); %counts up to half of AB block
            end
        elseif ss==3% last
            if  AB_border>=trgrange
            trgid=AB_border-trgrange+1:AB_border;
            else
                trgid=1:AB_border;
            end
        end
        
        if length(trgid)>size(tmptrg,1)
            trg=tmptrg(:,2);
        else
            trg=tmptrg(trgid,2);%sec
        end
        % --------------------------------------------------
        %%% cal rstar
        rst=ana_raster_ss(trg,win,spike);
        
        %%% cal psth
        psth = ana_peth_ss(rst,win,bin);
        
        %%% cal smooth psth
        bin_gs=100;%ms
        psth_gs=smoothing_gaussian_ss(psth,100,bin_gs,2);
        
        %%% cal z scored psth: version 2 JL 230414
        psth_gs_z=cal_z_score_v2(psth_gs,bin_baseline);
        % --------------------------------------------------
        history.AB_Rater{s,ss}=rst;
        history.AB_PSTH{s,ss}=psth(:);
        history.AB_PSTH_gs{s,ss}=psth_gs(:);
        history.AB_PSTH_gs_z{s,ss}=psth_gs_z(:);
        % --------------------------------------------------
        %%% detection of significant responses
        tmpS1=history.AB_PSTH_gs_z{s,ss};
        tmpS2=history.AB_PSTH_gs{s,ss};
        tmpS3=history.AB_PSTH_gs_z{s,ss};

        time_tmp=bin_odor;
        tmpd=tmpS1(time_tmp)>cri_z;
        [~,sig_id]=cal_continuous_num(tmpd,cri_bin1);  % z score more than 1.96 (p<0.05) or 2.58(p<0.01) for more than three consecutive bins (150ms) 
        tmp_sig1(s,ss)=sig_id;
        [history.AB_lat_peak_odor(s,ss),history.AB_lat_cb_odor(s,ss),history.AB_lat_30_odor(s,ss)]=cal_latency(tmpS1,time_tmp,cri_z,cri_bin1,cri_amp,bin_size);
        [history.AB_peak_odor(s,ss),history.AB_peak_pm_bins_odor(s,ss)]=cal_peak(tmpS2,time_tmp,bin_size,bin_range);     
        [history.AB_peak_odor_z(s,ss),history.AB_peak_pm_bins_odor_z(s,ss)]=cal_peak(tmpS3,time_tmp,bin_size,bin_range);     

        time_tmp=bin_delay;
        tmpd=tmpS1(time_tmp)>cri_z;
        [~,sig_id]=cal_continuous_num(tmpd,cri_bin1);
        tmp_sig2(s,ss)=sig_id;
        [history.AB_lat_peak_delay(s,ss),history.AB_lat_cb_delay(s,ss),history.AB_lat_30_delay(s,ss)]=cal_latency(tmpS1,time_tmp,cri_z,cri_bin2,cri_amp,bin_size);
        [history.AB_peak_delay(s,ss),history.AB_peak_pm_bins_delay(s,ss)]=cal_peak(tmpS2,time_tmp,bin_size,bin_range);     
        [history.AB_peak_delay_z(s,ss),history.AB_peak_pm_bins_delay_z(s,ss)]=cal_peak(tmpS3,time_tmp,bin_size,bin_range);     

        time_tmp=bin_rwd;
        tmpd=tmpS1(time_tmp)>cri_z;
        [~,sig_id]=cal_continuous_num(tmpd,cri_bin1);
        tmp_sig3(s,ss)=sig_id;
        [history.AB_lat_peak_rwd(s,ss),history.AB_lat_cb_rwd(s,ss),history.AB_lat_30_rwd(s,ss)]=cal_latency(tmpS1,time_tmp,cri_z,cri_bin1,cri_amp,bin_size);
        [history.AB_peak_rwd(s,ss),history.AB_peak_pm_bins_rwd(s,ss)]=cal_peak(tmpS2,time_tmp,bin_size,bin_range);     
        [history.AB_peak_rwd_z(s,ss),history.AB_peak_pm_bins_rwd_z(s,ss)]=cal_peak(tmpS3,time_tmp,bin_size,bin_range); 
        
        %%Add 230414
        time_tmp=bin_odor2;
        tmpd=tmpS1(time_tmp)>cri_z;
        [~,sig_id]=cal_continuous_num(tmpd,cri_bin1);
        tmp_sig4(s,ss)=sig_id;
        [history.AB_lat_peak_odor05_15(s,ss),history.AB_lat_cb_odor05_15(s,ss),history.AB_lat_30_odor05_15(s,ss)]=cal_latency(tmpS1,time_tmp,cri_z,cri_bin1,cri_amp,bin_size);
        [history.AB_peak_odor05_15(s,ss),history.AB_peak_pm_bins_odor05_15(s,ss)]=cal_peak(tmpS2,time_tmp,bin_size,bin_range);     
        [history.AB_peak_odor05_15_z(s,ss),history.AB_peak_pm_bins_odor05_15_z(s,ss)]=cal_peak(tmpS3,time_tmp,bin_size,bin_range);         
        
        time_tmp=bin_odor05_2;
        tmpd=tmpS1(time_tmp)>cri_z;
        [~,sig_id]=cal_continuous_num(tmpd,cri_bin1);
        tmp_sig5(s,ss)=sig_id;
        [history.AB_lat_peak_odor05_2(s,ss),history.AB_lat_cb_odor05_2(s,ss),history.AB_lat_30_odor05_2(s,ss)]=cal_latency(tmpS1,time_tmp,cri_z,cri_bin1,cri_amp,bin_size);
        [history.AB_peak_odor05_2(s,ss),history.AB_peak_pm_bins_odor05_2(s,ss)]=cal_peak(tmpS2,time_tmp,bin_size,bin_range);     
        [history.AB_peak_odor05_2_z(s,ss),history.AB_peak_pm_bins_odor05_2_z(s,ss)]=cal_peak(tmpS3,time_tmp,bin_size,bin_range); 
        
        time_tmp=bin_delay2;
        tmpd=tmpS1(time_tmp)>cri_z;
        [~,sig_id]=cal_continuous_num(tmpd,cri_bin1);
        tmp_sig6(s,ss)=sig_id;
        [history.AB_lat_peak_delay15_3(s,ss),history.AB_lat_cb_delay15_3(s,ss),history.AB_lat_30_delay15_3(s,ss)]=cal_latency(tmpS1,time_tmp,cri_z,cri_bin1,cri_amp,bin_size);
        [history.AB_peak_delay15_3(s,ss),history.AB_peak_pm_bins_delay15_3(s,ss)]=cal_peak(tmpS2,time_tmp,bin_size,bin_range);     
        [history.AB_peak_delay15_3_z(s,ss),history.AB_peak_pm_bins_delay15_3_z(s,ss)]=cal_peak(tmpS3,time_tmp,bin_size,bin_range); 
        %Added 230415
        time_tmp=bin_odor05_3;
        tmpd=tmpS1(time_tmp)>cri_z;
        [~,sig_id]=cal_continuous_num(tmpd,cri_bin1);
        tmp_sig7(s,ss)=sig_id;
        [history.AB_lat_peak_odor05_3(s,ss),history.AB_lat_cb_odor05_3(s,ss),history.AB_lat_30_odor05_3(s,ss)]=cal_latency(tmpS1,time_tmp,cri_z,cri_bin1,cri_amp,bin_size);
        [history.AB_peak_odor05_3(s,ss),history.AB_peak_pm_bins_odor03_2(s,ss)]=cal_peak(tmpS2,time_tmp,bin_size,bin_range);     
        [history.AB_peak_odor05_3_z(s,ss),history.AB_peak_pm_bins_odor05_3_z(s,ss)]=cal_peak(tmpS3,time_tmp,bin_size,bin_range); 
        
        for n=1:size(trgid,2) %Added 0.5-2s and 1.5-3s 230414
        
         rst1=rst(rst(:,2)==n,1); %picks up one trial
         psth_single = ana_peth_ss_single(rst1,win,0.05);
         psth_single=smoothing_gaussian_ss(psth_single,100,50,2);
        numSpike_pre1(n)=  mean(psth_single(bin_baseline));
        numSpike_odor1(n)= mean(psth_single(bin_odor2)); %Note- this odor period is 0.5-1.5s
        numSpike_delay(n)= mean(psth_single(bin_delay2)); %Note- this delay is 1.5-3s
        numSpike_rwd(n)= mean(psth_single(bin_rwd));
        
        numSpike_odor05_2(n) = mean(psth_single(bin_odor05_2)); %Added 230414
        
        numSpike_odor05_3(n) = mean(psth_single(bin_odor05_3)); %Added 230415
        end
        
        
        if ~isempty(numSpike_pre1) 
            [~,h1(s,ss)] = signrank(numSpike_pre1 ,numSpike_odor1);
            [h1_ttest(s,ss), ~] = ttest(numSpike_pre1 ,numSpike_odor1);

            [~,h5(s,ss)] = signrank(numSpike_pre1 ,numSpike_delay);
            
            [~,h6(s,ss)] = signrank(numSpike_pre1 ,numSpike_rwd);
            
            %Added 230414
            [~,h7(s,ss)] = signrank(numSpike_pre1 ,numSpike_odor05_2);
            %Added 230415
            [~,h8(s,ss)] = signrank(numSpike_pre1 ,numSpike_odor05_3);
            
        else
            h1(s,ss) = 0;
            h1_ttest(s,ss) = 0;
            h5(s,ss) = 0;
            h6(s,ss) = 0;
            h7(s,ss) = 0; 
            h8(s,ss) = 0;
        end

       
        rst1=[];
        numSpike_pre1=[];
        numSpike_odor1=[];
        numSpike_delay=[];
        numSpike_rwd=[];
        
        numSpike_odor05_2=[];
        numSpike_odor05_3=[];
        % --------------------------------------------------
    end
end


history.AB_SigID_odor_all =tmp_sig1(:,1); %For this original significance, odor = 0-1s
history.AB_SigID_delay_all=tmp_sig2(:,1); %Delay = 1-3s
history.AB_SigID_rwd_all  =tmp_sig3(:,1);
history.AB_SigID_odor_F   =tmp_sig1(:,2);
history.AB_SigID_delay_F  =tmp_sig2(:,2);
history.AB_SigID_rwd_F    =tmp_sig3(:,2);
history.AB_SigID_odor_L   =tmp_sig1(:,3);
history.AB_SigID_delay_L  =tmp_sig2(:,3);
history.AB_SigID_rwd_L    =tmp_sig3(:,3);

%Add 230414
history.AB_SigID_odor05_15_all =tmp_sig4(:,1); %Odor 0.5-1.5s
history.AB_SigID_odor05_2_all  =tmp_sig5(:,1); %Odor 0.5-2s
history.AB_SigID_delay15_3_all =tmp_sig6(:,1); %Delay 1.5-3s
history.AB_SigID_odor05_15_F   =tmp_sig4(:,2);
history.AB_SigID_odor05_2_F    =tmp_sig5(:,2);
history.AB_SigID_delay15_3_F   =tmp_sig6(:,2);
history.AB_SigID_odor05_15_L   =tmp_sig4(:,3);
history.AB_SigID_odor05_2_L    =tmp_sig5(:,3);
history.AB_SigID_delay15_3_L   =tmp_sig6(:,3);
%Add 230415
history.AB_SigID_odor05_3_all    =tmp_sig7(:,1);
history.AB_SigID_odor05_3_F    =tmp_sig7(:,2);
history.AB_SigID_odor05_3_L    =tmp_sig7(:,3);

history.AB_signrankOdor_all = h1(:,1); %NOTE: Signrank Odor period is bin_odor2 = 0.5-1.5s !
history.AB_signrankOdor_F = h1(:,2);
history.AB_signrankOdor_L = h1(:,3);

history.AB_signrankDelay_all = h5(:,1); %NOTE: Signrank Delay is bin_delay2 = 1.5-3s !
history.AB_signrankDelay_F = h5(:,2);
history.AB_signrankDelay_L = h5(:,3);

history.AB_signrankRwd_all = h6(:,1);
history.AB_signrankRwd_F = h6(:,2);
history.AB_signrankRwd_L = h6(:,3);

history.AB_ttestOdor_all = h1_ttest(:,1);
history.AB_ttestOdor_F = h1_ttest(:,2);
history.AB_ttestOdor_L = h1_ttest(:,3);

%Added 230414
history.AB_signrankOdor05_2_all = h7(:,1); %0.5-2s
history.AB_signrankOdor05_2_F = h7(:,2);
history.AB_signrankOdor05_2_L = h7(:,3);
%Added 230415
history.AB_signrankOdor05_3_all = h8(:,1); %0.5-2s
history.AB_signrankOdor05_3_F = h8(:,2);
history.AB_signrankOdor05_3_L = h8(:,3);

%%% cal diff - UNCHANGED 230414
clear tmp_sig1 tmp_sig2 tmp_sig3 tmp_sig4 tmp_sig5 tmp_sig6 tmp_sig7
for s=1:2
    history.AB_PSTH_gs_diff{s,1}=history.AB_PSTH_gs{s,2}-history.AB_PSTH_gs{s,3};
    % --------------------------------------------------
    %%% cal z scored psth
    psth_gs_z=cal_z_score_v2(history.AB_PSTH_gs_diff{s,1},bin_baseline);
    
    % --------------------------------------------------
    history.AB_PSTH_gs_diff_z{s,1}=psth_gs_z;

    % --------------------------------------------------
    %%% detection of significant responses
    tmpS1=history.AB_PSTH_gs_diff{s,1};

    time_tmp=bin_odor;
    tmpd=tmpS1(time_tmp)>cri_z;
    [~,sig_id]=cal_continuous_num(tmpd,cri_bin1);
    tmp_sig1(s,2)=sig_id;
    [history.AB_lat_peak_odor_diff(s,ss),history.AB_lat_cb_odor_diff(s,ss),history.AB_lat_30_odor_diff(s,ss)]=cal_latency(tmpS1,time_tmp,cri_z,cri_bin1,cri_amp,bin_size);
    [history.AB_peak_odor_diff(s,ss),history.AB_peak_pm_bins_odor_diff(s,ss)]=cal_peak(tmpS1,time_tmp,bin_size,bin_range);     

    time_tmp=bin_delay;
    tmpd=tmpS1(time_tmp)>cri_z;
    [~,sig_id]=cal_continuous_num(tmpd,cri_bin2);
    tmp_sig2(s,2)=sig_id;
    [history.AB_lat_peak_delay_diff(s,ss),history.AB_lat_cb_delay_diff(s,ss),history.AB_lat_30_delay_diff(s,ss)]=cal_latency(tmpS1,time_tmp,cri_z,cri_bin2,cri_amp,bin_size);
    [history.AB_peak_delay_diff(s,ss),history.AB_peak_pm_bins_delay_diff(s,ss)]=cal_peak(tmpS1,time_tmp,bin_size,bin_range);     

    time_tmp=bin_rwd;
    tmpd=tmpS1(time_tmp)>cri_z;
    [~,sig_id]=cal_continuous_num(tmpd,cri_bin1);
    tmp_sig3(s,3)=sig_id;
    [history.AB_lat_peak_rwd_diff(s,ss),history.AB_lat_cb_rwd_diff(s,ss),history.AB_lat_30_rwd_diff(s,ss)]=cal_latency(tmpS1,time_tmp,cri_z,cri_bin1,cri_amp,bin_size);
    [history.AB_peak_rwd_diff(s,ss),history.AB_peak_pm_bins_rwd_diff(s,ss)]=cal_peak(tmpS1,time_tmp,bin_size,bin_range);     
     % --------------------------------------------------
end
history.AB_SigID_odor_diff =tmp_sig1(:,1);
history.AB_SigID_delay_diff=tmp_sig2(:,2);
history.AB_SigID_rwd_diff  =tmp_sig3(:,3);





%% [3] General analysis - PSTH after novel odor presentation (ABCD block)
%%% config
trgrange=10;
% trial numbers for comparison between first and last X trials

%%% cal raster and psth
history.ABCD_name{1,1}='A all';
history.ABCD_name{1,2}=['A First',num2str(trgrange)];
history.ABCD_name{1,3}=['A Last',num2str(trgrange)];

history.ABCD_name{2,1}='B all';
history.ABCD_name{2,2}=['B First',num2str(trgrange)];
history.ABCD_name{2,3}=['B Last',num2str(trgrange)];

history.ABCD_name{3,1}='C all';
history.ABCD_name{3,2}=['C First',num2str(trgrange)];
history.ABCD_name{3,3}=['C Last',num2str(trgrange)];

history.ABCD_name{4,1}='D all';
history.ABCD_name{4,2}=['D First',num2str(trgrange)];
history.ABCD_name{4,3}=['D Last',num2str(trgrange)];

history.ABCD_name{1,4}=['A middle',num2str(trgrange)];
history.ABCD_name{2,4}=['B middle',num2str(trgrange)];
history.ABCD_name{3,4}=['C middle',num2str(trgrange)];
history.ABCD_name{4,4}=['D middle',num2str(trgrange)];

clear tmp_sig1 tmp_sig2 tmp_sig3
for s=1:4 % odor A to D
    oid=s; %odor id
    
    %%% get target time stamp
    trgid=history.time_stamp(:,3)==oid; %history.time_stamp created in Make_TimeStamp_list_ss; (:,3) points to odor ID list from Nlx data
    tmptrg=history.time_stamp(trgid,:);%sec
    tmptrg=tmptrg(~isnan(tmptrg(:,7)),:);% removing immature trials
    tmptrg=tmptrg(tmptrg(:,8)==0,:);% removing manual rwd trials
    
    if s==1 || s==2
        %AB_border=find(diff(tmptrg(:,4))<0); %original
        AB_border=find(tmptrg(:,1)>0,1,'first')-1;
        tmptrg=tmptrg(AB_border+1:end,:);
    end
    
    for ss=1:7% all, First, Last, middle,  Last20, 1/4, 3/4
        if ss==1% all
            trgid=1:size(tmptrg,1);
        elseif ss==2% first10
            trgid=1:trgrange;
        elseif ss==3% last10
            trgid=size(tmptrg,1)-trgrange+1:size(tmptrg,1);
        elseif ss==4% middle10
            trgid=round(size(tmptrg,1)/2)-trgrange/2+1:round(size(tmptrg,1)/2)+trgrange/2;
            if ~isempty(find(trgid<=0,1,'last'))  
                trgid=trgid(find(trgid<=0,1,'last')+1):trgid(find(trgid<=0,1,'last')+trgrange);     
            end
             %190724 m310 did not have enough AB trials for ABCD session
            if   length(trgid)>size(tmptrg,1)
                    disp('Not enough AB trials for ABCD session- reformatting trgid in cal_history_data line 320');
                    trgid=trgid(1:size(tmptrg,1));
            end
        elseif ss==5% last20
            if size(tmptrg,1)-trgrange*2+1>0
                trgid=size(tmptrg,1)-trgrange*2+1:size(tmptrg,1);
            else
                trgid=1:size(tmptrg,1);
            end
        elseif ss==6 % 1/3
            trgid=round(size(tmptrg,1)/3)-trgrange/2+1:round(size(tmptrg,1)/3)+trgrange/2;
        elseif ss==7 % 3/3
            
            trgid=round(size(tmptrg,1)*2/3)-trgrange/2+1:round(size(tmptrg,1)*2/3)+trgrange/2;
            
                
        end
        trg=tmptrg(trgid,2);%sec
        
        % --------------------------------------------------
        %%% cal rstar
        rst=ana_raster_ss(trg,win,spike);
        
        %%% cal psth
        psth = ana_peth_ss(rst,win,bin);
        
        %%% cal somooth psth
        bin_gs=50;%ms
        psth_gs=smoothing_gaussian_ss(psth,100,bin_gs,2);
        
        %%% cal z scored psth
        psth_gs_z=cal_z_score_v2(psth_gs,bin_baseline);
        % --------------------------------------------------
        history.ABCD_Rater{s,ss}=rst;
        history.ABCD_PSTH{s,ss}=psth;
        history.ABCD_PSTH_gs{s,ss}=psth_gs;
        history.ABCD_PSTH_gs_z{s,ss}=psth_gs_z;
        
        %         subplot(411);plot(rst(:,1),rst(:,2),'.')
        %         subplot(412);bar(psth)
        %         subplot(413);plot(psth_gs)
        %         subplot(414);plot(psth_gs_z)
        %         input('Enter')
        
        % --------------------------------------------------
        %%% detection of significant responses
        tmpS1=history.ABCD_PSTH_gs_z{s,ss};
        tmpS2=history.ABCD_PSTH_gs{s,ss};
        tmpS3=history.ABCD_PSTH_gs_z{s,ss};

        time_tmp=bin_odor;
        tmpd=tmpS1(time_tmp)>cri_z;
        [~,sig_id]=cal_continuous_num(tmpd,cri_bin1);
        tmp_sig1(s,ss)=sig_id;
        [history.ABCD_lat_peak_odor(s,ss),history.ABCD_lat_cb_odor(s,ss),history.ABCD_lat_30_odor(s,ss)]=cal_latency(tmpS1,time_tmp,cri_z,cri_bin1,cri_amp,bin_size);
        [history.ABCD_peak_odor(s,ss),history.ABCD_peak_pm_bins_odor(s,ss)]=cal_peak(tmpS2,time_tmp,bin_size,bin_range);     
        [history.ABCD_peak_odor_z(s,ss),history.ABCD_peak_pm_bins_odor_z(s,ss)]=cal_peak(tmpS3,time_tmp,bin_size,bin_range);     

        time_tmp=bin_delay;
        tmpd=tmpS1(time_tmp)>cri_z;
        [~,sig_id]=cal_continuous_num(tmpd,cri_bin2);
        tmp_sig2(s,ss)=sig_id;
        [history.ABCD_lat_peak_delay(s,ss),history.ABCD_lat_cb_delay(s,ss),history.ABCD_lat_30_delay(s,ss)]=cal_latency(tmpS1,time_tmp,cri_z,cri_bin2,cri_amp,bin_size);
        [history.ABCD_peak_delay(s,ss),history.ABCD_peak_pm_bins_delay(s,ss)]=cal_peak(tmpS2,time_tmp,bin_size,bin_range);     
        [history.ABCD_peak_delay_z(s,ss),history.ABCD_peak_pm_bins_delay_z(s,ss)]=cal_peak(tmpS3,time_tmp,bin_size,bin_range);     

        time_tmp=bin_rwd;
        tmpd=tmpS1(time_tmp)>cri_z;
        [~,sig_id]=cal_continuous_num(tmpd,cri_bin1);
        tmp_sig3(s,ss)=sig_id;
        [history.ABCD_lat_peak_rwd(s,ss),history.ABCD_lat_cb_rwd(s,ss),history.ABCD_lat_30_rwd(s,ss)]=cal_latency(tmpS1,time_tmp,cri_z,cri_bin1,cri_amp,bin_size);
        [history.ABCD_peak_rwd(s,ss),history.ABCD_peak_pm_bins_rwd(s,ss)]=cal_peak(tmpS2,time_tmp,bin_size,bin_range);     
        [history.ABCD_peak_rwd_z(s,ss),history.ABCD_peak_pm_bins_rwd_z(s,ss)]=cal_peak(tmpS3,time_tmp,bin_size,bin_range);
        
        %Added 230414
        time_tmp=bin_odor2;
        tmpd=tmpS1(time_tmp)>cri_z;
        [~,sig_id]=cal_continuous_num(tmpd,cri_bin1);
        tmp_sig4(s,ss)=sig_id;
        [history.ABCD_lat_peak_odor05_15(s,ss),history.ABCD_lat_cb_odor05_15(s,ss),history.ABCD_lat_30_odor05_15(s,ss)]=cal_latency(tmpS1,time_tmp,cri_z,cri_bin1,cri_amp,bin_size);
        [history.ABCD_peak_odor05_15(s,ss),history.ABCD_peak_pm_bins_odor05_15(s,ss)]=cal_peak(tmpS2,time_tmp,bin_size,bin_range);     
        [history.ABCD_peak_odor05_15_z(s,ss),history.ABCD_peak_pm_bins_odor05_15_z(s,ss)]=cal_peak(tmpS3,time_tmp,bin_size,bin_range);             
        
        time_tmp=bin_odor05_2;
        tmpd=tmpS1(time_tmp)>cri_z;
        [~,sig_id]=cal_continuous_num(tmpd,cri_bin1);
        tmp_sig5(s,ss)=sig_id;
        [history.ABCD_lat_peak_odor05_2(s,ss),history.ABCD_lat_cb_odor05_2(s,ss),history.ABCD_lat_30_odor05_2(s,ss)]=cal_latency(tmpS1,time_tmp,cri_z,cri_bin1,cri_amp,bin_size);
        [history.ABCD_peak_odor05_2(s,ss),history.ABCD_peak_pm_bins_odor05_2(s,ss)]=cal_peak(tmpS2,time_tmp,bin_size,bin_range);     
        [history.ABCD_peak_odor05_2_z(s,ss),history.ABCD_peak_pm_bins_odor05_2_z(s,ss)]=cal_peak(tmpS3,time_tmp,bin_size,bin_range);      
        
        time_tmp=bin_delay2;
        tmpd=tmpS1(time_tmp)>cri_z;
        [~,sig_id]=cal_continuous_num(tmpd,cri_bin2);
        tmp_sig6(s,ss)=sig_id;
        [history.ABCD_lat_peak_delay15_3(s,ss),history.ABCD_lat_cb_delay15_3(s,ss),history.ABCD_lat_30_delay15_3(s,ss)]=cal_latency(tmpS1,time_tmp,cri_z,cri_bin2,cri_amp,bin_size);
        [history.ABCD_peak_delay15_3(s,ss),history.ABCD_peak_pm_bins_delay15_3(s,ss)]=cal_peak(tmpS2,time_tmp,bin_size,bin_range);     
        [history.ABCD_peak_delay15_3_z(s,ss),history.ABCD_peak_pm_bins_delay15_3_z(s,ss)]=cal_peak(tmpS3,time_tmp,bin_size,bin_range);         
        
        %Add 230415
        time_tmp=bin_odor05_3;
        tmpd=tmpS1(time_tmp)>cri_z;
        [~,sig_id]=cal_continuous_num(tmpd,cri_bin1);
        tmp_sig7(s,ss)=sig_id;
        [history.ABCD_lat_peak_odor05_3(s,ss),history.ABCD_lat_cb_odor05_3(s,ss),history.ABCD_lat_30_odor05_3(s,ss)]=cal_latency(tmpS1,time_tmp,cri_z,cri_bin1,cri_amp,bin_size);
        [history.ABCD_peak_odor05_3(s,ss),history.ABCD_peak_pm_bins_odor05_3(s,ss)]=cal_peak(tmpS2,time_tmp,bin_size,bin_range);     
        [history.ABCD_peak_odor05_3_z(s,ss),history.ABCD_peak_pm_bins_odor05_3_z(s,ss)]=cal_peak(tmpS3,time_tmp,bin_size,bin_range);    
        
        for n=1:size(trgid,2)
        
        rst2=rst(rst(:,2)==n,1); %picks up one trial
         psth_single = ana_peth_ss_single(rst2,win,0.05);
         psth_single=smoothing_gaussian_ss(psth_single,100,50,2);
        numSpike_pre(n)=  mean(psth_single(bin_baseline));
        numSpike_odor(n)= mean(psth_single(bin_odor2));
        numSpike_delay(n)=mean(psth_single(bin_delay2));
        numSpike_rwd(n)=  mean(psth_single(bin_rwd));
        
        numSpike_odor05_2(n) = mean(psth_single(bin_odor05_2)); %Added 230414
        numSpike_odor05_3(n) = mean(psth_single(bin_odor05_3)); %Added 230415
        
        end
  
            [~,h2(s,ss)] = signrank(numSpike_pre ,numSpike_odor);
            [h2_ttest(s,ss), ~] = ttest(numSpike_pre ,numSpike_odor);
 
            [~,h3(s,ss)] = signrank(numSpike_pre ,numSpike_delay);
        
            [~,h4(s,ss)] = signrank(numSpike_pre ,numSpike_rwd);

            
            [~,h9(s,ss)] = signrank(numSpike_pre ,numSpike_odor05_2);
            [~,h10(s,ss)] = signrank(numSpike_pre ,numSpike_odor05_3);
            
        rst2=[];
        numSpike_pre=[];
        numSpike_odor=[];
        numSpike_delay=[];
        numSpike_rwd=[];
        numSpike_odor05_2=[];
        numSpike_odor05_3=[];
        % --------------------------------------------------
    end
end

history.ABCD_SigID_odor_all =tmp_sig1(:,1);
history.ABCD_SigID_delay_all=tmp_sig2(:,1);
history.ABCD_SigID_rwd_all  =tmp_sig3(:,1);
history.ABCD_SigID_odor_F   =tmp_sig1(:,2);
history.ABCD_SigID_delay_F  =tmp_sig2(:,2);
history.ABCD_SigID_rwd_F    =tmp_sig3(:,2);
history.ABCD_SigID_odor_L   =tmp_sig1(:,3);
history.ABCD_SigID_delay_L  =tmp_sig2(:,3);
history.ABCD_SigID_rwd_L    =tmp_sig3(:,3);

%Add 230414
history.ABCD_SigID_odor05_15_all =tmp_sig4(:,1);
history.ABCD_SigID_odor05_2_all=tmp_sig5(:,1);
history.ABCD_SigID_delay15_3_all  =tmp_sig6(:,1);
history.ABCD_SigID_odor05_15_F   =tmp_sig4(:,2);
history.ABCD_SigID_odor05_2_F  =tmp_sig5(:,2);
history.ABCD_SigID_delay15_3_F    =tmp_sig6(:,2);
history.ABCD_SigID_odor05_15_L   =tmp_sig4(:,3);
history.ABCD_SigID_odor05_2_L  =tmp_sig5(:,3);
history.ABCD_SigID_delay15_3_L    =tmp_sig6(:,3);
%Add 230415
history.ABCD_SigID_odor05_3_all  =tmp_sig7(:,1);
history.ABCD_SigID_odor05_3_F  =tmp_sig7(:,2);
history.ABCD_SigID_odor05_3_L  =tmp_sig7(:,3);


history.ABCD_signrankOdor_F = h2(:,2); %odor 0.5-1.5s
history.ABCD_signrankOdor_L = h2(:,3);
history.ABCD_signrankOdor_M = h2(:,4);
history.ABCD_signrankOdor_L20 = h2(:,5);

history.ABCD_signrankDelay_F = h3(:,2); %Delay 1.5-3s
history.ABCD_signrankDelay_L = h3(:,3);
history.ABCD_signrankDelay_M = h3(:,4);
history.ABCD_signrankDelay_L20 = h3(:,5);

history.ABCD_signrankRwd_F = h4(:,2);
history.ABCD_signrankRwd_L = h4(:,3);
history.ABCD_signrankRwd_M = h4(:,4);
history.ABCD_signrankRwd_L20 = h4(:,5);

history.ABCD_ttestOdor_F = h2_ttest(:,2);
history.ABCD_ttestOdor_L = h2_ttest(:,3);
history.ABCD_ttestOdor_M = h2_ttest(:,4);
history.ABCD_ttestOdor_L20 = h2_ttest(:,5);

%Add 230414
history.ABCD_signrankOdor05_2_F = h9(:,2); %odor 0.5-2s
history.ABCD_signrankOdor05_2_L = h9(:,3);
history.ABCD_signrankOdor05_2_M = h9(:,4);
history.ABCD_signrankOdor05_2_L20 = h9(:,5); 

%Add 230415
history.ABCD_signrankOdor05_3_F = h10(:,2); %odor 0.5-3s
history.ABCD_signrankOdor05_3_L = h10(:,3);
history.ABCD_signrankOdor05_3_M = h10(:,4);
history.ABCD_signrankOdor05_3_L20 = h10(:,5);


%% %%%%%%%%%%%%%%%%%%%%%  saving all cue ID and raster (for later analysis and shuffling - 200818 Kei)

% ABCD trials starting exactly at the first novel odor
history.odorIndex =[];
history.rst_all =[];
history.time_stamp2 = [];

     history.time_stamp2 = history.time_stamp(history.time_stamp(:,1)>0,:); % ABCD trials starting exactly at the first novel odor
     history.time_stamp2 = history.time_stamp2(~isnan(history.time_stamp2(:,7)),:); % removing immature trials
     history.time_stamp2 = history.time_stamp2(history.time_stamp2(:,8)==0,:);% removing manual rwd trials
     history.odorIndex   = history.time_stamp2(:,3); %%%%% odor ID for all trials
    
     
     trg2=history.time_stamp2(:,2);%sec
     history.rst_all=ana_raster_ss(trg2,win,spike); %%%%% raster for all trials (-2 to 6 sec)


     
% AB and ABCD trials bordered at LabVIEW scripts  
history.odorIndex_AB =[];
history.odorIndex_ABCD =[];

history.rst_all_AB =[];
history.rst_all_ABCD =[];

history.time_stamp2_AB = [];
history.time_stamp2_ABCD  = [];


starts = find(history.time_stamp(:,4)==1); % column 4 contains AB and ABCD trial IDs
startAB = starts(1);
startABCD = starts(2);



     history.time_stamp2_AB = history.time_stamp(startAB:startABCD-1,:); 
     history.time_stamp2_AB = history.time_stamp2_AB(~isnan(history.time_stamp2_AB(:,7)),:); % removing immature trials
     history.time_stamp2_AB = history.time_stamp2_AB(history.time_stamp2_AB(:,8)==0,:);% removing manual rwd trials
     history.odorIndex_AB   = history.time_stamp2_AB(:,3); %%%%% odor ID for all trials

     trg2_AB=history.time_stamp2_AB(:,2);%sec
     history.rst_all_AB=ana_raster_ss(trg2_AB,win,spike); %%%%% raster for all trials (-2 to 6 sec)


     history.time_stamp2_ABCD = history.time_stamp(startABCD:end,:); 
     history.time_stamp2_ABCD = history.time_stamp2_ABCD(~isnan(history.time_stamp2_ABCD(:,7)),:); % removing immature trials
     history.time_stamp2_ABCD = history.time_stamp2_ABCD(history.time_stamp2_ABCD(:,8)==0,:);% removing manual rwd trials
     history.odorIndex_ABCD   = history.time_stamp2_ABCD(:,3); %%%%% odor ID for all trials

     trg2_ABCD=history.time_stamp2_ABCD(:,2);%sec
     history.rst_all_ABCD=ana_raster_ss(trg2_ABCD,win,spike); %%%%% raster for all trials (-2 to 6 sec)





 %%%%%%%%%%%%%%%%%%%% 


%% cal_diff UNCHANGED
clear tmp_sig1 tmp_sig2 tmp_sig3 tmp_sig4 tmp_sig5 tmp_sig6
for s=1:4 %Odor A-D
    history.ABCD_PSTH_gs_diff{s,1}=history.ABCD_PSTH_gs{s,2}-history.ABCD_PSTH_gs{s,3};
    % --------------------------------------------------
    %%% cal z scored psth
    psth_gs_z=cal_z_score_v2(history.ABCD_PSTH_gs_diff{s,1},bin_baseline);
    % --------------------------------------------------
    history.ABCD_PSTH_gs_diff_z{s,1}=psth_gs_z;
    
    %     subplot(211);plot(history.ABCD_PSTH_gs_diff{s,1})
    %     subplot(212);plot(psth_gs_z)
    %     input('Enter')
    
    % --------------------------------------------------
    %%% detection of significant responses
    tmpS1=history.ABCD_PSTH_gs_diff{s,1};


    time_tmp=bin_odor;
    tmpd=tmpS1(time_tmp)>cri_z;
    [~,sig_id]=cal_continuous_num(tmpd,cri_bin1);
    tmp_sig1(s,1)=sig_id;
    [history.ABCD_lat_peak_odor_diff(s,ss),history.ABCD_lat_cb_odor_diff(s,ss),history.ABCD_lat_30_odor_diff(s,ss)]=cal_latency(tmpS1,time_tmp,cri_z,cri_bin1,cri_amp,bin_size);
    [history.ABCD_peak_odor_diff(s,ss),history.ABCD_peak_pm_bins_odor_diff(s,ss)]=cal_peak(tmpS1,time_tmp,bin_size,bin_range);     

    time_tmp=bin_delay;
    tmpd=tmpS1(time_tmp)>cri_z;
    [~,sig_id]=cal_continuous_num(tmpd,cri_bin2);
    tmp_sig2(s,2)=sig_id;
    [history.ABCD_lat_peak_delay_diff(s,ss),history.ABCD_lat_cb_delay_diff(s,ss),history.ABCD_lat_30_delay_diff(s,ss)]=cal_latency(tmpS1,time_tmp,cri_z,cri_bin2,cri_amp,bin_size);
    [history.ABCD_peak_delay_diff(s,ss),history.ABCD_peak_pm_bins_delay_diff(s,ss)]=cal_peak(tmpS1,time_tmp,bin_size,bin_range);     

    time_tmp=bin_rwd;
    tmpd=tmpS1(time_tmp)>cri_z;
    [~,sig_id]=cal_continuous_num(tmpd,cri_bin1);
    tmp_sig3(s,3)=sig_id;
    [history.ABCD_lat_peak_rwd_diff(s,ss),history.ABCD_lat_cb_rwd_diff(s,ss),history.ABCD_lat_30_rwd_diff(s,ss)]=cal_latency(tmpS1,time_tmp,cri_z,cri_bin1,cri_amp,bin_size);
    [history.ABCD_peak_rwd_diff(s,ss),history.ABCD_peak_pm_bins_rwd_diff(s,ss)]=cal_peak(tmpS1,time_tmp,bin_size,bin_range);     
    % --------------------------------------------------
end
history.ABCD_SigID_odor_diff =tmp_sig1(:,1);
history.ABCD_SigID_delay_diff=tmp_sig2(:,2);
history.ABCD_SigID_rwd_diff  =tmp_sig3(:,3);









