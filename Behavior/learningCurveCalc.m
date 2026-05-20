% consider learning curve at only matureful trials
% add the part to calculate the mature rate around each matureful trial
% as well
% I prefer to calculate the LC or matureRate by placing the trial in the
% center, such as i-10:i+10. But lab exisiting code use trials after
% current idx. I will follow what lab does
% Li YUAN, 2026-Mar
function [LC,matureRate] = learningCurveCalc(trialInfo,LC_window_all,LC_window_each)

matureNum = sum(trialInfo.mature == 1);
matureIdx = find(trialInfo.mature == 1);
matureOdor = trialInfo.odor(trialInfo.mature == 1);
matureOutcome = trialInfo.outcome(trialInfo.mature == 1);
if any(isnan(matureOutcome))
    error('Trial outcome error')
end

%% calculate based on all odors
if matureNum > LC_window_all
    LC.ABCD = nan(6,matureNum-LC_window_all+1);
    matureRate.all = nan(1,matureNum-LC_window_all+1);
    for i = 1:(matureNum-LC_window_all)

        Go_1 = sum(matureOdor(i:i+LC_window_all-1) == 1); % # of odor A presentations in this 20 trials
        NoGo_1 = sum(matureOdor(i:i+LC_window_all-1) == 2); % # of odor B... etc...
        Go_2 = sum(matureOdor(i:i+LC_window_all-1) == 3);
        NoGo_2 = sum(matureOdor(i:i+LC_window_all-1) == 4);

        Hit_1 = sum((matureOdor(i:i+LC_window_all-1) == 1) & (matureOutcome(i:i+LC_window_all-1) == 1));
        CR_1  = sum((matureOdor(i:i+LC_window_all-1) == 2) & (matureOutcome(i:i+LC_window_all-1) == 4));
        Hit_2 = sum((matureOdor(i:i+LC_window_all-1) == 3) & (matureOutcome(i:i+LC_window_all-1) == 1));
        CR_2  = sum((matureOdor(i:i+LC_window_all-1) == 4) & (matureOutcome(i:i+LC_window_all-1) == 4));

        %LC stands for Learning Curve
        LC.ABCD(1,i)=((Hit_1 + CR_1)/(Go_1 + NoGo_1)) * 100; % Correct rate AB
        LC.ABCD(2,i)= Hit_1/Go_1 * 100; % Hit rate for A
        LC.ABCD(3,i)= CR_1/NoGo_1 * 100; % CR rate for B
        LC.ABCD(4,i)=((Hit_2 + CR_2)/(Go_2 + NoGo_2)) * 100; % Correct rate CD
        LC.ABCD(5,i)= Hit_2/Go_2 * 100; % Hit rate for C
        LC.ABCD(6,i)= CR_2/NoGo_2 * 100; % CR rate for D


        matureEnd = min([matureIdx(i)+LC_window_all,trialInfo.num]);
        matureLabel = trialInfo.mature(matureIdx(i):matureEnd);
        matureRate.all(i) = sum(matureLabel)/length(matureLabel)*100;
    end

    % fill last 20 trials with the same data
    lastWindow = (matureNum - LC_window_all + 1);
    %disp
    Go_1 = sum(matureOdor(lastWindow:matureNum) == 1); % # of odor A presentations in this 20 trials
    NoGo_1 = sum(matureOdor(lastWindow:matureNum) == 2); % # of odor B... etc...
    Go_2 = sum(matureOdor(lastWindow:matureNum) == 3);
    NoGo_2 = sum(matureOdor(lastWindow:matureNum) == 4);

    Hit_1 = sum((matureOdor(lastWindow:matureNum) == 1) & (matureOutcome(lastWindow:matureNum) == 1));
    CR_1  = sum((matureOdor(lastWindow:matureNum) == 2) & (matureOutcome(lastWindow:matureNum) == 4));
    Hit_2 = sum((matureOdor(lastWindow:matureNum) == 3) & (matureOutcome(lastWindow:matureNum) == 1));
    CR_2  = sum((matureOdor(lastWindow:matureNum) == 4) & (matureOutcome(lastWindow:matureNum) == 4));

    %LC stands for Learning Curve
    LC.ABCD(1,lastWindow:matureNum)=((Hit_1 + CR_1)/(Go_1 + NoGo_1)) * 100; % Correct rate AB
    LC.ABCD(2,lastWindow:matureNum)= Hit_1/Go_1 * 100; % Hit rate for A
    LC.ABCD(3,lastWindow:matureNum)= CR_1/NoGo_1 * 100; % CR rate for B
    LC.ABCD(4,lastWindow:matureNum)=((Hit_2 + CR_2)/(Go_2 + NoGo_2)) * 100; % Correct rate CD
    LC.ABCD(5,lastWindow:matureNum)= Hit_2/Go_2 * 100; % Hit rate for C
    LC.ABCD(6,lastWindow:matureNum)= CR_2/NoGo_2 * 100; % CR rate for D
    
    matureEnd = min([matureIdx(matureNum - LC_window_all + 1)+LC_window_all,trialInfo.num]);
    matureLabel = trialInfo.mature(matureIdx(i):matureEnd);
    matureRate.all(matureNum - LC_window_all + 1) = sum(matureLabel)/length(matureLabel)*100;

else
    fprintf('Less than %d matureful trials, no LC calculation\n',LC_window_all);
    LC.ABCD = nan(6,matureNum);
    matureRate.all = nan(1,matureNum);
end

%% calculate based on each odors
mature_A_outcome = matureOutcome(matureOdor == 1);
mature_B_outcome = matureOutcome(matureOdor == 2);
mature_C_outcome = matureOutcome(matureOdor == 3);
mature_D_outcome = matureOutcome(matureOdor == 4);

if length(mature_A_outcome) >= LC_window_each
    LC.A = smooth_movWin(mature_A_outcome==1, LC_window_each)*100;
else
    LC.A = nan;
end

if length(mature_B_outcome) >= LC_window_each
    LC.B = smooth_movWin(mature_B_outcome==4, LC_window_each)*100;
else
    LC.B = nan;
end

if length(mature_C_outcome) >= LC_window_each
    LC.C = smooth_movWin(mature_C_outcome==1, LC_window_each)*100;
else
    LC.C = nan;
end

if length(mature_D_outcome) >= LC_window_each
    LC.D = smooth_movWin(mature_D_outcome==4, LC_window_each)*100;
else
    LC.D = nan;
end

%% calculate AB and CD
mature_AB_outcome = matureOutcome(matureOdor==1 | matureOdor==2);
mature_CD_outcome = matureOutcome(matureOdor==3 | matureOdor==4);
AB_correct = zeros(size(mature_AB_outcome));
AB_correct(mature_AB_outcome==1 | mature_AB_outcome==4) = 1;
if length(mature_AB_outcome) >= LC_window_each
    LC.AB = smooth_movWin(AB_correct, LC_window_each)*100;
else
    LC.AB = nan;
end

CD_correct = zeros(size(mature_CD_outcome));
CD_correct(mature_CD_outcome==1 | mature_CD_outcome==4) = 1;
if length(mature_CD_outcome) >= LC_window_each
    LC.CD = smooth_movWin(CD_correct, LC_window_each)*100;
else
    LC.CD = nan;
end

end