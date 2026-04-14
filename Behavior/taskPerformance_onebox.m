% this code analyze the go-nogo task with cue
% LiYUAN, 2026-Mar
function trial = taskPerformance_onebox(data,Fs,setting,pairId)
% premature trial = 0, successful trial = 1;
% odor 1A, 2B, 3C, 4D and so on % Or odd number is go, even number is nogo
% trial success or not, result in gngGetTS_check_error_AB_ABCD_labView: 0 fail, 1 success
% in successful trials, evnt.action in gngGetTS_check_error_AB_ABCD_labView: (1: Hit, 2: Miss, 3: FalseAlarm, 4: CorrectRejection)
% 1 & 2 for A,C, 3 & 4 for B,D
% correct, 1&4, error, 2&3
PLOT = 1;
% voltage threshold for each line to be considered on
odorThres = 2;
LEDThres = 1;
lickThres = 1.5;
rewardThres = 2;
allodors = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R']; % make sure go is always odd and nogo is always even

odorLenMax = 5; % unit sec
delayAddup = 2; % unit sec, in case lab view give a bit more time

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
    h = figure;
    h.Position = [100,100,900,900];
    f(1) = subplot(4,1,1);
    ts = (1:size(data,2))/Fs;
    plot(ts,data(setting.odorA_ch,:))
    hold on
    plot(ts,data(setting.odorB_ch,:))
    plot(ts,data(setting.odorC_ch,:))
    plot(ts,data(setting.odorD_ch,:))
    % plot(ts,data(setting.odorE_ch,:))
    % plot(ts,data(setting.odorF_ch,:))
    plot([0,ts(end)],[odorThres,odorThres],'k')
    text(0,odorThres+0.5,'Thres')
    legend({'A','B','C','D'})
    ylabel('Odor Amp (V)')

    f(2) = subplot(4,1,2);
    plot(ts,data(setting.led_ch,:))
    hold on
    plot([0,ts(end)],[LEDThres,LEDThres],'k')
    text(0,LEDThres+0.5,'Thres')
    ylabel('LED Amp (V)')
    ylim([0 4])

    f(3) = subplot(4,1,3);
    plot(ts,data(setting.lick_ch,:))
    hold on
    plot([0,ts(end)],[lickThres,lickThres],'k')
    text(0,lickThres+0.5,'Thres')
    ylabel('Lick Amp (V)')
    title('Lick signal not used in deciding correct or not')
    ylim([0 5])

    f(4) = subplot(4,1,4);
    plot(ts,data(setting.sucrose_ch,:))
    hold on
    plot(ts,data(setting.quinine_ch,:))
    plot([0,ts(end)],[rewardThres,rewardThres],'k')
    text(0,rewardThres+0.5,'Thres')
    legend({'S','Q'})
    ylabel('Reward Amp (V)')
    xlabel('Time (s)')

    linkaxes(f,'x')

end

odor_sig.A = data(setting.odorA_ch,:) > odorThres;
odor_sig.B = data(setting.odorB_ch,:) > odorThres;
odor_sig.C = data(setting.odorC_ch,:) > odorThres; % same channel used for C,E,G ETC
odor_sig.D = data(setting.odorD_ch,:) > odorThres; % same channel used for D,F,H ETC
% odor_sig.E = data(setting.odorE_ch,:) > 2;
% odor_sig.F = data(setting.odorF_ch,:) > 2;
sucrose_sig = data(setting.sucrose_ch,:) > rewardThres;
quinine_sig = data(setting.quinine_ch,:) > rewardThres;
led_sig = data(setting.led_ch,:) > LEDThres;
lick_sig = data(setting.lick_ch,:) > lickThres;
delay_length = setting.delay_len;
response_length = setting.response_len;


% get total trial Num

trial.num = 0;
onsetIndTemp = [];
offsetIndTemp = [];
% odorTemp = []; % 1A,2B,3C,4D,5E,6F etc
odorTemp = []; % 1A,2B,3 novel (C,E,G ETC), 4 novel (D,F,H ETC)

% detect odor onset and odor name
for k = 1:length(odorName)

    thisOdor = upper(odorName{k});
    [~,odorID] = ismember(thisOdor,allodors);
    if odorID == 0
        error('Odor name detection is wrong')
    else
        if odorID > 2 % A B is always 1, 2
            if rem(odorID,2) == 0
                thisOdor = 'C';
                odorID = 3;
            else
                thisOdor = 'D';
                odorID = 4;
            end
        end

        odorPulse = odor_sig.(thisOdor);
        pulseDiff = diff(odorPulse);
        pulseDiff = [0;pulseDiff(:)];
        trialNumTemp = sum(pulseDiff==1);
        trial.num = trial.num + trialNumTemp;
        onsetIndTemp = [onsetIndTemp;find(pulseDiff==1)];
        offsetIndTemp = [offsetIndTemp;find(pulseDiff==-1)];
        odorTemp = [odorTemp;odorID * ones(trialNumTemp,1)];
    end
end

% sort by time
[trial.onsetInd,idx] = sort(onsetIndTemp);
trial.odor = odorTemp(idx);

% correct trial num
% remove last trial if last onset is too close to end which is shorter than
% delay + response time
% maxLim = ceil(trial.Onset(end) + Fs*(delay_length + response_length));
% if maxLim > length(trial.Onset)
%     trial.Num = trial.Num - 1;
%     trial.Onset = trial.Onset(1:end-1);
%     trial.Odor = trial.Odor(1:end-1);
% end

trial.offsetInd = nan(trial.num,1);
trial.success = nan(trial.num,1);
trial.outcome = nan(trial.num,1);

% get the outcome of each trial
for j = 1:trial.num-1
    % whether a trial is successful
    % check whether there is LED after lick onset
    % or to check whether there is licking during delay

    % find the offset
    offsetInd = offsetIndTemp(offsetIndTemp > trial.onsetInd(j));
    trial.offsetInd(j) = min(offsetInd);
    if trial.offsetInd(j) - trial.onsetInd(j) > ceil(odorLenMax*Fs) % assume no odor lasts for more than 5 sec
        error('Odor duration error')
    end

    sigWindow = trial.offsetInd(j) : (trial.offsetInd(j) + ceil(Fs*(delay_length + response_length + delayAddup)));
    if sum(led_sig(sigWindow)) == 0  % LED not on, failed
        trial.success(j) = 0;
    else
        trial.success(j) = 1;
        % decide the outcome
        if rem(trial.odor(j),2) == 1  % go trials, odd number in labeling
            if sum(sucrose_sig(sigWindow)) > 0
                trial.outcome(j) = 1;
            else
                trial.outcome(j) = 2;
            end
        elseif rem(trial.odor(j),2) == 0  % no go trials, even number in labeling
            if sum(quinine_sig(sigWindow)) > 0
                trial.outcome(j) = 3;
            else
                trial.outcome(j) = 4;
            end
        else
            % pause
            error('Odor is wrong')
        end

    end
end

% take the last trial consider boundry issues
% 1. if the odor onset is too close to end
j = trial.num;
offsetInd = offsetIndTemp(offsetIndTemp > trial.onsetInd(j));
trial.offsetInd(j) = min(offsetInd);
if trial.offsetInd(j) - trial.onsetInd(j) > odorLenMax*Fs % assume no odor lasts for more than 5 sec
    error('Odor duration error')
end

maxLen = size(data,2);
trialLim = trial.offsetInd(j) + ceil(Fs*(delay_length + response_length + 0.2));
if maxLen >= trialLim
    sigWindow = trial.offsetInd(j) : (trial.offsetInd(j) + ceil(Fs*(delay_length + response_length + 0.2)));
    if sum(led_sig(sigWindow)) == 0  % LED not on, failed
        trial.success(j) = 0;
    else
        trial.success(j) = 1;
        % decide the outcome
        if trial.odor(j) == 1 || trial.odor(j) == 3  % go trials
            if sum(sucrose_sig(sigWindow)) > 0
                trial.outcome(j) = 1;
            else
                trial.outcome(j) = 2;
            end
        elseif trial.odor(j) == 2 || trial.odor(j) == 4  % no go trials
            if sum(quinine_sig(sigWindow)) > 0
                trial.outcome(j) = 3;
            else
                trial.outcome(j) = 4;
            end
        else
            % pause
            error('Odor is wrong')
        end
    end
else
    trial.num = trial.num-1;
    trial.onsetInd = trial.onsetInd(1:trial.num);
    trial.offsetInd = trial.offsetInd(1:trial.num);
    trial.odor = trial.odor(1:trial.num);
    trial.success = trial.success(1:trial.num);
    trial.outcome = trial.outcome(1:trial.num);
end