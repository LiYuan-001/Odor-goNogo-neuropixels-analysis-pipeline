% =========================================================
% Return counts of each nidq channel type that compose
% the timepoints stored in binary file.
%
function [XA,DW,SY] = ChannelCountsOBX(meta)
    M = str2num(meta.snsXaDwSy);
    XA = M(1);
    DW = M(2);
    SY = M(3);
end % ChannelCountsNI