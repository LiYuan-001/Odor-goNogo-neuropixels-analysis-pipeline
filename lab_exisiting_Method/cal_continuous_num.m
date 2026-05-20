function [continous_bin,id]=cal_continuous_num(data,c_num)

%data=[0 1 0 0 1 0 0 1 0 1 1 0 0 1 1 1 0 1 1 0 1 1 1 1 1];
%c_num=4;

tmp=data;
for i=1:length(tmp)-c_num+1
    x(i)=sum(tmp(i:c_num-1+i));
end

continous_bin=find(x>=c_num);
if isempty(continous_bin)
    continous_bin=NaN;
end

if sum(continous_bin)>1
    id=1;
else
    id=0;
end
end



% subplot(221)
% plot(tmp,'o-');hold on
% % plot([0 length(tmp)+1],[1 1],'r:')
% ylim([-0.5 1.5])
%
%  subplot(222)
%  plot(cumsum(tmp))
%
%  subplot(223)
% for i=1:length(tmp)-c_num+1
%    x(i)=sum(tmp(i:c_num-1+i));
% end
%  plot(x)
%