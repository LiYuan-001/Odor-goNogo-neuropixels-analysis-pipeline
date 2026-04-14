% =========================================================
% Return gain for ith channel stored in the obx file.
%
% ichan is a saved channel index, rather than an original
% (acquired) index.
%
function gain = ChanGainOBX(ichan, sa,  meta)
    if ichan <= sa
        gain = str2double(meta.niMNGain);
    else
        gain = 1;
    end
end % ChanGainNI