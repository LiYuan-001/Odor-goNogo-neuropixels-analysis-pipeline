function dataArray = GainCorrectNI(dataArray, chanList, meta)

    [MN,MA] = ChannelCountsNI(meta);
    fI2V = Int2Volts(meta);

    for i = 1:length(chanList)
        j = chanList(i);    % index into timepoint
        conv = fI2V / ChanGainNI(j, MN, MA, meta);
        dataArray(j,:) = dataArray(j,:) * conv;
    end
end % GainCorrectNI

function fI2V = Int2Volts(meta)
    if strcmp(meta.typeThis, 'imec')
        if isfield(meta,'imMaxInt')
            maxInt = str2num(meta.imMaxInt);
        else
            maxInt = 512;
        end
        fI2V = str2double(meta.imAiRangeMax) / maxInt;
    elseif strcmp(meta.typeThis, 'nidq')
        fI2V = str2double(meta.niAiRangeMax) / 32768;
    elseif strcmp(meta.typeThis, 'obx')
        fI2V = str2double(meta.obAiRangeMax) / 32768;
    end
end % Int2Volts