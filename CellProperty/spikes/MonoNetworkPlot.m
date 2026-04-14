function MonoNetworkPlot(sig_con,cellID,cellLabel)

weight = zeros(length(cellID),length(cellID));
allcel = unique(sig_con(:));
cellLabel2 = zeros(length(allcel),1);
for i = 1:size(sig_con,1)
    weight(sig_con(i,1),sig_con(i,2)) = 1;
    weight(sig_con(i,2),sig_con(i,1)) = 1;
end
counts = sum(weight~=0);
counts2 = counts*2+1;
G = graph(weight);
LWidths = 2*G.Edges.Weight/max(G.Edges.Weight);
h = plot(G,'NodeLabel',cellID,'LineWidth',LWidths);
for nn = 1:length(cellID)
    if cellLabel(nn) == 1
        highlight(h,nn, 'MarkerSize', counts2(nn),'NodeColor','r','Marker','^')
    elseif cellLabel(nn) == 2
        highlight(h,nn, 'MarkerSize', counts2(nn),'NodeColor','b','Marker','o')
    end
end
end
