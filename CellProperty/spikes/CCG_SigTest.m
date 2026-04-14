function [Pval,Pred,Bounds,sig_con,sig_con_inh] = CCG_SigTest(ccgR,tR,conv_w,sigWindow,binSize,alpha)
% get CI for each CCG
Pval=nan(length(tR),2,2);
Pred=zeros(length(tR),2,2);
Bounds=zeros(size(ccgR,1),2,2);
sig_con = [0,0,0];
sig_con_inh = [0,0,0];
%     sig_con1 = [];
%     TruePositive = nan(cellNum,cellNum);
%     FalsePositive = nan(cellNum,cellNum);
Pcausal = nan(2,2);

cch = ccgR(:,1,2);
centerbins = ceil(length(cch)/2);
%             cch(centerbins) = NaN;
[pvals,pred,qvals]=bz_cch_conv(cch,conv_w);

% Store predicted values and pvalues for subsequent plotting
Pred(:,1,2)=pred;
Pval(:,2,1)=pvals(:);
Pred(:,1,2)=flipud(pred(:));
Pval(:,2,1)=flipud(pvals(:));

nBonf = round(sigWindow/binSize)*2;

hiBound=poissinv(1-alpha/nBonf,pred);
loBound=poissinv(alpha/nBonf, pred);
Bounds(:,1,2,1)=hiBound;
Bounds(:,1,2,2)=loBound;
Bounds(:,2,1,1)=flipud(hiBound(:));
Bounds(:,2,1,2)=flipud(loBound(:));

% sig = cch>hiBound | cch < loBound;
sig = cch>hiBound;

% Find if significant periods falls in monosynaptic window +/- 4ms
prebins = floor(length(cch)/2 - sigWindow/binSize):floor(length(cch)/2);
postbins = floor(length(cch)/2)+1:floor(length(cch)/2 + sigWindow/binSize);
cchud  = flipud(cch);
sigud  = flipud(sig);
% sigpost=max(cch(postbins))>poissinv(1-alpha,max(cch(prebins)));
% sigpre=max(cchud(postbins))>poissinv(1-alpha,max(cchud(prebins)));

%define likelihood of being a connection
pvals_causal = 1 - poisscdf( max(cch(postbins)) - 1, max(cch(prebins) )) - poisspdf( max(cch(postbins)), max(cch(prebins)  )) * 0.5;
pvals_causalud = 1 - poisscdf( max(cchud(postbins)) - 1, max(cchud(prebins) )) - poisspdf( max(cchud(postbins)), max(cchud(prebins)  )) * 0.5;

%can go negative for very small p-val - beyond comp. sig. dig

if pvals_causalud<0
    pvals_causalud = 0;
end

if pvals_causal<0
    pvals_causal = 0;
end

Pcausal(1,2) = pvals_causal;
Pcausal(2,1) = pvals_causalud;

% %check which is bigger
% if (any(sigud(prebins)) && sigpre)
%     %test if causal is bigger than anti causal
%     sig_con = [1,1,2];
% end
% if any(sig(postbins)) && sigpost
%     sig_con = [1,2,1];
% end
if any(sig(prebins))
    sig_con = [1,2,1];
end
if any(sig(postbins))
    sig_con = [1,1,2];
end
    
% % % % % % % % % % % % % % % % % % % % % % %
% INHIBITORY
sig_inh = cch<loBound;
cchud = flipud(cch);
sigud_inh = flipud(sig_inh);
% sigpost_inh = min(cch(postbins))<poissinv(alpha,max(cch(prebins)));
% sigpre_inh = min(cchud(postbins))<poissinv(alpha,max(cchud(prebins)));

% % check which is bigger
% if (any(sigud_inh(postbins)) && sigpre_inh)
%     %test if causal is bigger than anti causal
%     sig_con_inh = [1,1,2];
% end
% 
% if any(sig_inh(postbins)) && sigpost_inh
%     sig_con_inh = [1,2,1];
% end

if any(sig_inh(prebins))
    sig_con_inh = [1,2,1];
end
if any(sig_inh(postbins))
    sig_con_inh = [1,1,2];
end
end