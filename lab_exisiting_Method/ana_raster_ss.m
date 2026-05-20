function rst = ana_raster_ss(trg,win,spike)

rst=[];
win1=win(1);
win2=win(2);
trg=trg(trg>=spike(1) & trg<=spike(end));
nodewndw=[trg+win1 trg+win2];
for k=1:size(nodewndw,1)
    count=sum(spike<=nodewndw(k,2))-sum(spike<=nodewndw(k,1));
    if count~=0
        for l=1:count
            t=intersect(find(spike>=nodewndw(k,1)),find(spike<=nodewndw(k,2)));%
            rst(end+1,1)=spike(t(l))-trg(k,1); %column1: when did the firing happen in respect to the start time of themovement(trigger‚©‚çspike‚Ü‚Å‚ÌŽžŠÔ(ms))
            rst(end,2)=k; %column2: which movement does the firing relate with(heavy or light,go or nogoŽŽ?s‚Ì’Ê‚µ”Ô?†?j           
        end
    else %there is no firing during the up movement being assessed
        rst(end+1,1)=NaN;
        rst(end,2)=k;
    end
end
if isempty(rst)
    rst=[0 0];
end
