% Having acquired a block of raw obx data using ReadBin(),
% convert values voltages. The conversion is only applied 
% to the saved-channel indices in chanList. The conversion
% factor is the same for all channels, because the gain is fixed.
% Remember saved-channel indices are in range [1:nSavedChans].
% The dimensions of the dataArray remain unchanged. ChanList
% examples:
%
%   [2,6,20]    % just these three channels
%
%
function dataArray = GainCorrectOBX(dataArray, chanList, meta)

    fI2V = Int2Volts(meta);

    for i = 1:length(chanList)
        j = chanList(i);    % index into timepoint
        dataArray(j,:) = dataArray(j,:) * fI2V;
    end
end % GainCorrectOBX