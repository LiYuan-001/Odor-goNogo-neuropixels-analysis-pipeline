% map sync2 time according to sync1
% use pulses match the sync frequency
% Li Yuan, Jul-23-2022, UCSD
% 
function timeStamps = syncTime2(sync_1,sync_2,sync1_Ts,Fs,pulseWidth)

Fs_1 = 1./(sync1_Ts(2)-sync1_Ts(1));
timeStamps = zeros(length(sync_2),1);
rising_sync_Fil = zeros(length(sync_1),1);
rising_sync2_Fil = zeros(length(sync_2),1);

sync_Diff = diff(sync_1);
if sum(sync_Diff>0)>0
    rising_sync_Ind = find(sync_Diff>0);
    rising_sync_Ind = rising_sync_Ind + 1;
else
    error('No sync signal.')
end


sync2_Diff = diff(sync_2);
if sum(sync2_Diff>0)>0
    rising_sync2_Ind = find(sync2_Diff>0);
    rising_sync2_Ind = rising_sync2_Ind + 1;
else
    error('No sync2 signal.')
end

% get only piulseWidth rising edges & remove short intervals
temp = find(round(diff(sync1_Ts(rising_sync_Ind)),1) ~= pulseWidth);
% temp = temp+1;
% rising_sync_Ind(temp) = [];
rising_sync_Ind(temp+1) = [];
i = 1;
while i < length(rising_sync_Ind)
    if round((sync1_Ts(rising_sync_Ind(i+1)) - sync1_Ts(rising_sync_Ind(i))),2) == pulseWidth
        rising_sync_Fil(rising_sync_Ind(i)) = 1;
        rising_sync_Fil(rising_sync_Ind(i+1)) = 1;
    end 
    i = i+1;
end

% rising_sync_Ind(temp) = [];
% rising_sync_Fil(rising_sync_Ind(~temp)) = 1;

sync2_Ts_Temp = (0:length(sync_2))./Fs;
temp = find(round(diff(sync2_Ts_Temp(rising_sync2_Ind)),1) ~= pulseWidth);
% temp = temp+1;
% rising_sync2_Ind(temp) = [];
rising_sync2_Ind(temp+1) = [];
i = 1;
while i < length(rising_sync2_Ind)
    if round((sync2_Ts_Temp(rising_sync2_Ind(i+1)) - sync2_Ts_Temp(rising_sync2_Ind(i))),2) == pulseWidth
        rising_sync2_Fil(rising_sync2_Ind(i)) = 1;
        rising_sync2_Fil(rising_sync2_Ind(i+1)) = 1; 
    end
    i = i+1;
end

% rising_sync2_Ind(temp) = [];
% rising_sync2_Fil(rising_sync2_Ind(~temp)) = 1;

clear sync_1 sync_2 sync2_Ts_Temp

% detect rising pulse based on filtered signals
sync_Diff = diff(rising_sync_Fil);
rising_sync_Ind = find(sync_Diff>0);
rising_sync_Ind = rising_sync_Ind + 1;

sync2_Diff = diff(rising_sync2_Fil);
rising_sync2_Ind = find(sync2_Diff>0);
rising_sync2_Ind = rising_sync2_Ind + 1;

% insert steps for some specific recordings that ttl is missed
if any(round(unique(diff(rising_sync_Ind))/pulseWidth/Fs_1)>1)
    step = rising_sync_Ind(2)-rising_sync_Ind(1);
    temp = diff(rising_sync_Ind);
    indTemp = find(temp>pulseWidth*1.2*Fs_1);
    for k = 1:length(indTemp)
        insertCount = round(temp(indTemp(k))/Fs_1/pulseWidth);
        for m = 1:insertCount-1
            rising_sync_Fil(rising_sync_Ind(indTemp(k))+step*m) = 1;
        end
    end
    sync_Diff = diff(rising_sync_Fil);
    rising_sync_Ind = find(sync_Diff>0);
    rising_sync_Ind = rising_sync_Ind + 1;
end

if any(round(unique(diff(rising_sync2_Ind))/pulseWidth/Fs)>1)
    step = rising_sync2_Ind(2)-rising_sync2_Ind(1);
    temp = diff(rising_sync2_Ind);
    indTemp = find(temp>pulseWidth*1.2*Fs);
    for k = 1:length(indTemp)
        insertCount = round(temp(indTemp(k))/Fs/pulseWidth);
        for m = 1:insertCount-1
            rising_sync2_Fil(rising_sync2_Ind(indTemp(k))+step*m) = 1;
        end
    end
    sync2_Diff = diff(rising_sync2_Fil);
    rising_sync2_Ind = find(sync2_Diff>0);
    rising_sync2_Ind = rising_sync2_Ind + 1;
end
    
if length(rising_sync_Ind) ~= length(rising_sync2_Ind)
    error('Pulse sizes are different')
end

if length(unique(round(diff(sync1_Ts(rising_sync_Ind)),2))) > 1
    error('missing sync ind')
end

% map timestamps for every sample in sync2 based on sample amounts in two
% rising edge

% condition 1
% from beginning to first rising edge
startInd = 1;
endInd = rising_sync2_Ind(1);
startTs = sync1_Ts(1);
endTs = sync1_Ts(rising_sync_Ind(1));
sampleSizeTemp = endInd - startInd;    
tsFill = (endTs - startTs)/sampleSizeTemp;
timeStamps(startInd:endInd) = startTs:tsFill:endTs;
% condition 2
% between two rising edges
for i = 1:length(rising_sync_Ind)-1
    startInd = rising_sync2_Ind(i);
    endInd = rising_sync2_Ind(i+1);
    
    startTs = sync1_Ts(rising_sync_Ind(i));
    endTs = sync1_Ts(rising_sync_Ind(i+1));
    
    sampleSizeTemp = endInd - startInd;
    if sampleSizeTemp > 1.5*Fs*pulseWidth || sampleSizeTemp < 0.5*Fs*pulseWidth
        error('Pulse detection was wrong')
    end
    
    tsFill = (endTs - startTs)/sampleSizeTemp;
    timeStamps(startInd:endInd) = startTs:tsFill:endTs;
end
% condition 3
% from last rising edge to the end
% keep using previous tsFill
startInd = rising_sync2_Ind(end);
endInd = length(rising_sync2_Fil);
startTs = sync1_Ts(rising_sync_Ind(end));
sampleSizeTemp = endInd - startInd;
timeStamps(startInd:endInd) = startTs:tsFill:(startTs+tsFill*sampleSizeTemp);  

clear sync_Diff rising_sync_Ind sync2_Diff rising_sync2_Ind rising_sync2_Fil rising_sync_Fil sync_1 sync_2

% double check
% diff ts should be same if it comes to round 3 decimals
tsDiff = diff(timeStamps);
decimalPart = floor(log10(1/Fs));
if length(unique(round(tsDiff,decimalPart))) > 1
    error('TimeSTamp assignment is wrong')
end