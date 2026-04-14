% =========================================================
% Return sample rate as double.
%
function srate = SampRate(meta)
    if strcmp(meta.typeThis, 'imec')
        srate = str2double(meta.imSampRate);
    elseif strcmp(meta.typeThis, 'nidq')
        srate = str2double(meta.niSampRate);
    elseif strcmp(meta.typeThis, 'obx')
        srate = str2double(meta.obSampRate);
    else
        error('Unrecognized type')
    end
end % SampRate