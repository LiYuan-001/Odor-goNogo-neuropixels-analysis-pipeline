
function [smoothed]=smooth_movWin(data, windowSize)
% smooth usign sliding window of 3 bins

if size(data,1) > size(data,2)
    data=data'; 
end

totalSize = size(data,2);
smoothed = nan(1,totalSize);

for i=1:(totalSize-windowSize+1) %average using moving window
   smoothed(i)= mean(data(i:i+windowSize-1));
   
end

for j=(totalSize-windowSize+2):totalSize    %process tips   
    smoothed(j)= mean(data(j:totalSize));
end

end
