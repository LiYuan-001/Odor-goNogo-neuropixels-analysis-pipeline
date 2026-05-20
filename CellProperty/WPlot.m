function WPlot(W,kth,plotWflat)
% plots seqNMF factors W
set(gcf, 'color', 'w');

if nargin<3
    plotWflat = 1;
end

[N,K,L] = size(W); 
color_palet = [[0 .6 .3]; [.7 0 .7]; [1 .6 0];  [.1 .3 .9];  [1 .1 .1];  [0 .9 .3]; [.4 .2 .7]; [.7 .2 .1]; [.1 .8 1 ]; [1 .3 .7]; [.2 .8 .2]; [.7 .4 1]; [.9 .6 .4]; [0 .6 1]; [1 .1 .3]]; 
color_palet = repmat(color_palet, ceil(K/size(color_palet,1)),1); 
kColors = color_palet(1:K,:); 
% cleanerW = W.*(W>max(W(:))*.1);
% cleanerH = H.*(H>max(H(:))*.1); 
%% set widths of subplots
m = .1; % margin
ww = .2; % width of W plot
if plotWflat
    wwflat = .05; % width of Wflat plot
else 
    wwflat = 0;
end

wdata = 1-ww-wwflat-2*m; 
sep = ceil(L*.5); 


%% plot W's
% axW = subplot('Position', [m m ww hdata]);cla
% hold on
set(gca, 'ColorOrder', kColors); 

WsToPlot = zeros(N,K*(L+sep)); 
XsToPlot = zeros(3,K); 
YsToPlot = [zeros(1,K); N*ones(2,K)]+.5;
for ki = 1:K
    XsToPlot(:,ki) = [(L+sep)*(ki-1) (L+sep)*(ki-1) (L+sep)*(ki-1)+L+1];
    WsToPlot(:,((L+sep)*(ki-1)+1):((L+sep)*(ki-1)+L)) = squeeze(W(:,ki,:));    
end
plot(XsToPlot,YsToPlot, 'linewidth', 2);
hold on
plot(XsToPlot(:,kth),YsToPlot(:,kth), 'k','linewidth', 4);
clims = [0 prctile(WsToPlot(WsToPlot>0),99)]; % if all W's are empty this line will bug
helper.grayscalePatchPlot(flipud(WsToPlot));%, clims(2)); 

xlim([0 K*(L+sep)]);

set(gca, 'ydir', 'normal')
axis off


% %% plot Wflat (collapse out L dimension of W)
% if plotWflat
%     axWflat = subplot('Position', [m+ww+wdata m wwflat hdata]);cla
%     hold on
%     set(gca, 'ColorOrder', kColors); 
%     plot(squeeze(sum(W,3)), N:-1:1,'>', 'markersize', 2.5);
%     axis tight
% 
%     xlims = xlim; 
%     xlim([xlims(2)*.1 xlims(2)])
%     set(gca, 'ydir', 'normal')
%     axis off
% end

