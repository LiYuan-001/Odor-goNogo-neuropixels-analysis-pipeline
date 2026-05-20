function psth = ana_peth_ss(rst,win,bin)

win1=win(1);
win2=win(2);

if sum(rst)==0
    psth=zeros(1,length(win1:bin:win2),1);
else
    if ~isempty(rst)
        trials=length(unique(rst(:,2)));
        psth=histc(rst(:,1),win1:bin:win2)/trials/bin;
        
        % here the unit is converted to count/second!
    else
        psth=zeros(1,length(win1:bin:win2));
    end
end
psth=psth(:);
end
