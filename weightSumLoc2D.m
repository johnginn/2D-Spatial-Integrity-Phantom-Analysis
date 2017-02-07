% Function to calculated the location of a circle based on the weighted sum
% of signals within a given image.
%
% Input:
% xSearch Array of x-locations used to calculate the lcoation of the cylinder
% ySearch Array of y-locations used to calculate the location of the cylinder
% sigImg The volume containing the signals (must match x-,y-,z-Search)
%
% Output:
% weightCoord = [weightXInd, weightYInd];
%
% John Ginn
% Created: 11/2/16
% Modified: 11/2/16

function [weightCoord] = ...
    weightSumLoc2D(xSearch,ySearch,sigImg)

% calculate the total signal
totSig = 0;
for xStep = 1:length(xSearch)
    for yStep = 1:length(ySearch)
            totSig = totSig + sigImg(yStep,xStep);
    end
end

% normalize the signal
normSig = sigImg./totSig;
weightXInd = 0;
weightYInd = 0;
for xStep = 1:length(xSearch)
    for yStep = 1:length(ySearch)
            currentSig = normSig(yStep,xStep);
            % weight the locations by the signal
            weightXInd = weightXInd + currentSig*xSearch(xStep);
            weightYInd = weightYInd + currentSig*ySearch(yStep);
    end
end

weightCoord = [weightXInd, weightYInd];

end