function ratio = burstiness(ts,p)

totalNum = length(ts);
spkInterval = diff(ts);
burstNum = sum(spkInterval <= p.burst);
ratio = burstNum/totalNum*100;

end