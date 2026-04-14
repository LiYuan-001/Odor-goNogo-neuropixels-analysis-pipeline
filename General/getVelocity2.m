function [v,vHat] = getVelocity2(x,y,t)

dy = diff(y);
dx = diff(x);
dt = diff(t);

v = arrayfun(@(u,v)norm([u v]),dy,dx)./dt;
if size(v,1) > 1
    v = [v(1);v];
else
    v = [v(1),v];
end
%mean velocity
vHat = mean(v);
end