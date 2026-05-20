function match = findContinuousOdorMatch(onebox, labview_event)
% findContinuousOdorMatch
% Finds the longest continuous matching odor sequence between:
%   trialInfo_pair2.odor where trialInfo_pair2.mature == 1
% and
%   evnt_pair2.odorID
%
% Output:
%   match.bestLen              length of matched sequence
%   match.oneboxStart          start index in oneboxOdor
%   match.oneboxEnd            end index in oneboxOdor
%   match.eventStart           start index in evnt_pair2.odorID
%   match.eventEnd             end index in evnt_pair2.odorID
%   match.matchedTrialInfoIdx  original indices in trialInfo_pair2
%   match.matchedEventIdx      indices in evnt_pair2.odorID
%   match.matchedOneboxOdor    matched odor sequence from trialInfo
%   match.matchedEventOdor     matched odor sequence from event

    oneBoxInd  = onebox.mature == 1;
    oneboxOdor = onebox.odor(oneBoxInd);
    eventOdor  = labview_event.odorID;

    oneboxOdor = oneboxOdor(:);
    eventOdor  = eventOdor(:);

    oneBoxSuccessIdx = find(oneBoxInd);

    bestLen = 0;
    bestOneboxStart = NaN;
    bestEventStart  = NaN;

    for i = 1:numel(oneboxOdor)
        for j = 1:numel(eventOdor)

            k = 0;

            while (i+k <= numel(oneboxOdor)) && ...
                  (j+k <= numel(eventOdor)) && ...
                  isequal(oneboxOdor(i+k), eventOdor(j+k))

                k = k + 1;
            end

            if k > bestLen
                bestLen = k;
                bestOneboxStart = i;
                bestEventStart  = j;
            end
        end
    end

    if bestLen == 0
        match = struct();
        match.bestLen = 0;
        match.oneboxStart = NaN;
        match.oneboxEnd = NaN;
        match.eventStart = NaN;
        match.eventEnd = NaN;
        match.matchedTrialInfoIdx = [];
        match.matchedEventIdx = [];
        match.matchedOneboxOdor = [];
        match.matchedEventOdor = [];
        warning('No continuous odor match found.');
        return
    end

    bestOneboxEnd = bestOneboxStart + bestLen - 1;
    bestEventEnd  = bestEventStart  + bestLen - 1;

    matchedTrialInfoIdx = oneBoxSuccessIdx(bestOneboxStart:bestOneboxEnd);
    matchedEventIdx     = bestEventStart:bestEventEnd;

    match = struct();
    match.bestLen = bestLen;

    match.oneboxStart = bestOneboxStart;
    match.oneboxEnd   = bestOneboxEnd;

    match.eventStart = bestEventStart;
    match.eventEnd   = bestEventEnd;

    match.matchedTrialInfoIdx = matchedTrialInfoIdx;
    match.matchedEventIdx     = matchedEventIdx;

    match.matchedOneboxOdor = oneboxOdor(bestOneboxStart:bestOneboxEnd);
    match.matchedEventOdor  = eventOdor(bestEventStart:bestEventEnd);
end