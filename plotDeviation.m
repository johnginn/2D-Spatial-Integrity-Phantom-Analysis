% Function to plot the deviation of the cylinders with respect to their
% distance from siocenter. The location of the errorbar is plotted at the
% center of the bin.
%
% Input:
% distance (mm) The distance the cylinder is from isocenter
% deviation (mm) The deviation of the cylinder location from ground truth
% binSize (mm) The bin size for grouping the data with respect to their distance from isocenter 
%
% Output:
% figHandle The handle to the figure
%
% John Ginn
% Created: 7/8/16
% Modified: 7/15/16

function [figHandle] = plotDeviation(distance,deviation,binSize)
maxDist = max(distance);
nBins = ceil(maxDist/binSize); % the number of bins
binPlot = (binSize.*(1:(nBins - 1)) - 0.5*binSize); % locate the bins at the middle of the bin
% make the last bin be located at the center of the bin, the last bin may
% be smaller than the rest if the max(distance) is not a multiple of the
% binSize -->(furthest datapoint - boundary of closest other bin)*0.5 +  boundary of closest other bin
% binPlot(nBins) = (max(distance) - binSize*(nBins-1))*0.5 + binSize*(nBins-1); % plot all bins at the center of the bin
binPlot(nBins) = (binSize*nBins - 0.5*binSize); % equally sized bins
binBounds = binSize.*(0:nBins); % the bounds for the bins
% sort the data into the different bins
cylindersInBin = zeros(1,nBins); % count the number of cylinders in each bin
currentBin = 1; % the index of the current bin
currentBinNum = zeros(length(distance),1);
% initialize the string array for the box plot by making an array of empty
% strings
maxStrLength = 0; % length of the longest string
for step = 1:length(binPlot)
    if maxStrLength < length(num2str(binPlot(step)))
        % new longest string
        maxStrLength = length(num2str(binPlot(step)));
    end
end
% make the empty string array
maxStr = blanks(maxStrLength);
for step = 1:length(distance)
   binLocation(step,:) = maxStr;  
end

for step = 1:length(distance)
    currentDist = distance(step);
    currentDev = deviation(step);
    % reset whether or not the bin has been found
    binFound = 0; % bin found? y = 1, n = 0
    binStep = 1; % reset cycling through the bins
    while binFound == 0;
        if ((currentDist >= binBounds(binStep))&&...
                (currentDist <= binBounds(binStep + 1)))
            % the bin has been found
            binFound = 1;
            % add one count to the bin
            cylindersInBin(binStep) = cylindersInBin(binStep) + 1;
            % add the cylinder and deviation to the bin
            bincylinders{cylindersInBin(binStep),binStep} = currentDist;
            binDeviation{cylindersInBin(binStep),binStep} = currentDev;
            % the current bin location, store for boxplot
            currentBinString = num2str(binPlot(binStep)); 
            currentStringLength = length(currentBinString);
            currentBinNum(step) = binPlot(binStep); % for sorting
            % avoiding string length mismatch...
            binLocation(step,:) = [currentBinString blanks(maxStrLength - currentStringLength)]; 
        else
            % cylinder not in current bin
            binStep = binStep + 1;
        end
    end
end
% sort the bins into the correct order for plotting the box-plot



% calculate the mean and standard deviation of the cylinder deviation from
% ground truth for each bin
totDev = 0;
binMeanDev = zeros(1,nBins);
binStdv = zeros(1,nBins);
for stepBin = 1:nBins
    ncylinders = cylindersInBin(stepBin); % number of cylinders in the bin
    deviationArray = zeros(1,ncylinders); % array of cylinder deviation
    totDev = 0; % reset the total deviation
    for stepcylinders = 1:ncylinders
        % store the individual deviation
        deviationArray(stepcylinders) = binDeviation{stepcylinders,stepBin};
        % the total deviatiion
        totDev = totDev + deviationArray(stepcylinders);
    end
    % calculate the mean deviation
    if ncylinders == 0;
       ncylinders = 1; % avoid dividing by zero 
    end
    binMeanDev(stepBin) = totDev/ncylinders;
    binStdv(stepBin) = std(deviationArray);
    
end

% sort the bins into the proper order
[sortedBins, sortIndex] = sort(currentBinNum);
for sortStep = 1:length(sortIndex)
   sortedDeviation(sortStep) = deviation(sortIndex(sortStep)); 
   sortedBinLocation(sortStep,:) = binLocation(sortIndex(sortStep),:);
end
whiskerSize = 1000; % make sure there are no outliers, specifies maximum whisker length

figure; % make a new figure so this one doesn't get replaced
h = boxplot(sortedDeviation,sortedBinLocation,'Whisker',whiskerSize);
lineWidth = 1.2;
set(h,{'linew'},{lineWidth})
xlabel('Distance from isocenter (mm)','FontSize',22)
ylabel('2D Deviation (mm)','FontSize',22)
title('Deviation from Ground Truth','FontSize',22)

% the old plot showing all the individual data
% % plot the data
% % black and white
% plot(distance,deviation,'.','MarkerFaceColor',0.7.*[1 1 1])
% hold on
% plot(binPlot,binMeanDev,'ok','MarkerSize',8)
% hold on
% errorbar(binPlot,binMeanDev,binStdv,'k','LineWidth',1.8);% color version
% % plot(distance,deviation,'.b')
% % hold on
% % errorbar(binPlot,binMeanDev,binStdv,'r','LineWidth',1.8);
% % hold on
% % plot(binPlot,binMeanDev,'or','MarkerSize',8)
% ylim([0 (max(deviation)*1.1)])
% xlim([0 max(distance)*1.05])
% xlabel('Distance from Isocenter (mm)','FontSize',22)
% ylabel('Deviation (mm)','FontSize',22)
% title('Deviation from Ground Truth','FontSize',22)
% legend('Data','Bin Average')


end