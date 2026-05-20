
function [smoothed2]=smooth_gaussian2(data,sigma,bin_width,dimension)
%sigma bin_width in sec (.004 for 4ms, bin_width=.001 for 1000Hz)
% rate2S=smooth_gaussian2(rate2,100,bin,1);
% dimension==1: smooth across columns for data in each rows
% dimension==2: smooth across raws for data in each columns

ratio=10*sigma/bin_width; 


if dimension==1
    
    for i=1:size(data,1)
        data2=[repmat(data(i,1),1,ratio), data(i,:), repmat(data(i,size(data,2)),1,ratio)]; % prolong edges
        edges=[-3*sigma:bin_width:3*sigma];
        kernel=normpdf(edges,0,sigma);
        kernel=kernel*bin_width; %multiply by bin width
        center=ceil(length(edges)/2);
        smoothed=conv(data2,kernel);
        smoothed=smoothed(center:numel(data2)+center-1);
        smoothed2(i,:)=smoothed(ratio+1:numel(data2)-ratio); % cut edges to proper length
        
    end
    
    
    
elseif  dimension==2
    
    for i=1:size(data,2)
        data2=[repmat(data(1,i),ratio,1); data(:,i); repmat(data(size(data,1),i),ratio,1)];
        edges=[-3*sigma:bin_width:3*sigma];
        kernel=normpdf(edges,0,sigma);
        kernel=kernel*bin_width; %multiply by bin width
        center=ceil(length(edges)/2);
        smoothed=conv(data2,kernel);
        smoothed=smoothed(center:numel(data2)+center-1);
        smoothed2(:,i)=smoothed(ratio+1:numel(data2)-ratio);
        
    end
    
    
end

end
