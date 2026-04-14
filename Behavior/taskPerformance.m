function trial = taskPerformance(data,Fs,setting,pairId)

% premature trial = 0, successful trial = 1;
% odor 1A, 2B, 3C, 4D, 5E, 6F
% trial success or not, result in gngGetTS_check_error_AB_ABCD_labView: 0 fail, 1 success
% in successful trials, evnt.action in gngGetTS_check_error_AB_ABCD_labView: (1: Hit, 2: Miss, 3: FalseAlarm, 4: CorrectRejection)
% 1 & 2 for A,C, 3 & 4 for B,D
% correct, 1&3, error, 2&4
PLOT = 1;
% voltage threshold for each line to be considered on
odorThres = 2;
LEDThres = 1;
lickThres = 2;
rewardThres = 2;
allodors = ['A','B','C','D','E','F'];


if contains(pairId,'pair1')
    odorName = setting.pair1;
elseif contains(pairId,'pair2')
    odorName = setting.pair2;
else
    error('Odor pair not decided')
end

odorNum = length(odorName);

% get signal and digitize

if PLOT == 1
    h = figure(1);
    subplot(3,1,1)
    h.Position = [100,100,600,900];
    ts = (1:size(data,2))/Fs;
    plot(ts,data(setting.odorA_ch,:))
    hold on
    plot(ts,data(setting.odorB_ch,:))
    plot(ts,data(setting.odorC_ch,:))
    plot(ts,data(setting.odorD_ch,:))
    % plot(ts,data(setting.odorE_ch,:))
    % plot(ts,data(setting.odorF_ch,:))
    plot(ts,odorThres*ones(size(data,2),1),'k')
    legend({'A','B','C','D','E','F'})
    ylabel('Odor Amp (V)')

    subplot(3,1,2)
    plot(ts,data(setting.LEDThres,:))
    ylabel('LED Amp (V)')

end

odor_sig.A = data(setting.odorA_ch,:) > 2;
odor_sig.B = data(setting.odorB_ch,:) > 2;
odor_sig.C = data(setting.odorC_ch,:) > 2;
odor_sig.D = data(setting.odorD_ch,:) > 2;
odor_sig.E = data(setting.odorE_ch,:) > 2;
odor_sig.F = data(setting.odorF_ch,:) > 2;
sucrose_sig = data(setting.sucrose_ch,:) > volt_high;
quinine_sig = data(setting.quinine_ch,:) > volt_high;
led_sig = data(setting.led_ch,:) > volt_high;
lick_sig = data(setting.lick_ch,:) > volt_high;
delay_length = setting.delay_len;
response_length = setting.response_len;


% get total trial Num

trial.Num = 0;
trial.Onset = [];
trial.Odor = []; % 1A,2B,3C,4D,5E,6F

% detect odor onset and odor name
for k = 1:length(odorName)

    thisOdor = upper(odorName{k});
    odorPulse = odor_sig.(thisOdor);
    odorPulse = odorPulse(:);
    pulseDiff = [0;odorPulse];
    trialNumTemp = sum(pulseDiff==1);
    trial.Num = trial.Num + trialNumTemp;
    trial.Onset = [trial.Onset;pulseDiff==1];


    [~,odorID] = ismember(odorName,allodors);
    if odorID == 0
        error('Odor name is wrong')
    else
        trial.Odor = [trial.Odor;odorID * ones(trialNumTemp,1)];
    end
end

% correct trial num
% remove last trial if last onset is too close to end which is shorter than
% delay + response time
% maxLim = ceil(trial.Onset(end) + Fs*(delay_length + response_length));
% if maxLim > length(trial.Onset)
%     trial.Num = trial.Num - 1;
%     trial.Onset = trial.Onset(1:end-1);
%     trial.Odor = trial.Odor(1:end-1);
% end

trial.success = nan(trial.Num,1);
trial.outcome = nan(trial.Num,1);
% get the outcome of each trial
for j = 1:trial.Num-1
    % whether a trial is successful
    % check whether there is LED after lick onset
    % or to check whether there is licking during delay
    if sum(led_sig(trial.Onset(j) + Fs*(delay_length + response_length))) == 0  % LED not on, failed
        trial.success(j) = 0;
    else
        trial.success(j) = 1;

        % decide the outcome

        
    end
end
end