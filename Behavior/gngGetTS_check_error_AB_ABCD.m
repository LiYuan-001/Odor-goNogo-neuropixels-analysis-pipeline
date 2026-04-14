%% gngGetTS_check_error_AB_ABCD

% extracts "event" that contains odorID and behavioral action in each trial

% 1) premature licking (0: premature, 1: success)
% 2) odor on time
% 3) odor off time
% 4) odor ID
% ---------------
% 5) Results (1: Hit, 2: Miss, 3: FalseAlarm, 4: CorrectRejection)
% 6) (Go)licking time  (NoGo)responce time end
% 7) Manual water (0: NoMW, 1<: MW)


%% check errors added by LabView in first few trials (based on shogo's experience)
fTrial=events(1:7);
if sum(fTrial(1:2)<10)==2% error type A (If first two events do not include a timestamp for some reason- added an extra starting value)
    n_start=2;disp('Error A: Immature or Mature') %Start reading events at entry #2 (assumes #1 was randomly added)   
    ErrID=1;
elseif sum(fTrial([1 7])<10)==2% no error, mature trial
    n_start=1;disp('No Error');ErrID=0;
elseif fTrial(1)<10 && fTrial(7)>10
    if fTrial(1)==0% no error, premature trial
        n_start=1;disp('No Error');ErrID=0;
    elseif fTrial(1)==1 %  error type 2
        n_start=6;disp('Error Type2')
        ErrID=1;  
    end
else
    Error=input('New Type Error');
end

t_type=[];
signal_ss=[];
% manual_rwd=[];
i = 1; %event
j = 1; %results
while n_start + 3 <= length(events)
    %n_start
    if events(n_start) == 0 % premature-licking trial
        results.suc_delay(j) = 0;
        results.mw(j) = NaN;
        t_type(j) = events(n_start + 3);
        results.action(j) = NaN;
        results.RT(j) = NaN;
         
        % [1]odor type ,[2]odor onset, [3]NaN, [4]NaN
        tmp_ss = [events(n_start + 3) events(n_start + 2) NaN NaN];
        signal_ss=[signal_ss; tmp_ss];
        

        j = j+1;
        n_start = n_start + 5;
    elseif events(n_start) == 1 % sucsess trial
      % evnt.odorON(i) = events(n_start + 2);   Shogo's code; should be +1
        evnt.odorON(i) = events(n_start + 1);       
        evnt.odorID(i) = events(n_start + 3);
        evnt.action(i) = events(n_start + 4); % 1:Hit, 2:Miss, 3:FalseAlearm, 4: Correct Rejection
        if evnt.action(i)==1 || evnt.action(i)==4
            CoEr=1;
        else
            CoEr=0;
        end
        if (evnt.odorID(i) == 1)&&(evnt.action(i)==3)
            evnt.action(i)=2;
        end

        evnt.RT(i) = events(n_start + 5) - events(n_start + 2);
        evnt.MW(i) = events(n_start + 6);
        
        results.suc_delay(j) = 1;
        results.mw(j) = evnt.MW(i);
        t_type(j) = evnt.odorID(i);
        results.RT(j) = evnt.RT(i);
 
        % [1]odor type ,[2]odor onset, [3]Correct / Miss, [4]manual_rwd
        tmp_ss = [events(n_start + 3) events(n_start + 2) CoEr evnt.MW(i);];
        signal_ss=[signal_ss; tmp_ss];

        if (evnt.action(i) == 1)||(evnt.action(i) == 3)
            results.action(j) = 1; % go trial
        else
            results.action(j) = 0; % nogo trial
        end
        
        if i > 1
            if (evnt.action(i-1) == 2) || (evnt.action(i-1) == 3)
                evnt.correction(i) = 1; % correction trial (prev trial was incorrect)
            else
                evnt.correction(i) = 0;
            end
        end
        
        n_start = n_start + 7;
        i = i+1;
        j = j+1;
    end     
end

