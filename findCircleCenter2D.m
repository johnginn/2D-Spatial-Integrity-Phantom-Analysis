% This function searches for the center of the cylinder in the phantom given
% the volume data of the phantom an array containing the volume of the 
% cylinder. 
%
% Input:
% imData The image of signal data (dimensions: x,y)
% circleData The image of the cylinder (dimensions: x,y)
% startLocation [xInd,yInd] The coordinates of the first (top left) cylinder
% voxelHeight (mm/voxel) The height of the pixel
% voxelWidth (mm/voxel) The width of the pixel
% searchDist (mm) distance in each direction to search for the cylinder using 
% the correlation coefficinet method
%
% searchDistWeight (mm) distance in each direction from the location found
% by the correlation coeffienct method used to calculate the weighted sum position
%
% upScale y = 1, n = 0 search for the cylinder center based off of upscaled data
% cylinderBoundary  y = 1, n = 0 Whether or not a boundary of 0's surrounds the cylinder model
% MRorCT ('mr' or 'ct') Determines whether the template contrast exhibits a CT or MR scan 
% upscaleFactor The factor the image is upscaled by
%
% Output:
% coordCylinderData [xInd,yInd] The coordinates of the cylininder as the appear to user in imshow
% 
% coordCylinderPlot [yInd,xInd] The coordinates of the cylinder 
% for plotting because the coordinates are reversed in imshow
%
% optCorrelation The optimal correlation coefficient
%
% weightCoord [xInd,yInd] The coordinates determined by the weighted sum method
% 
% weightCoordPlot [yInd,xInd] The coordinates determined by the weighted sum method
%
% coordcylinderDataUp [xInd,yInd] The coordinates of the cylinder from upscaled
% analysis as the appear to user in imshow
% 
% coordcylinderPlotUp [yInd,xInd] Nearst voxel in original image using
% the upscaled cylinder location
%
% optCorrelationUp The optimal correlation coefficient from upscaled
% analysis
%
% normX The x-indicies of the volume searched
% normY The y-indicies of the volume searched
%
% John Ginn
% Created: 11/2/16
% Modified: 11/4/16

function [coordCylinderData,coordCylinderPlot, optCorrelation,...
    weightCoord, weightCoordPlot,coordCylinderDataUp, coordCylinderPlotUp,...
    optCorrelationUp, normX, normY] = ...
    findCircleCenter2D(imgData, circleData,startLocation,voxelHeight,voxelWidth,...
    searchDist,searchDistWeight,upScale,cylinderBoundary,MRorCT,upscaleFactor)
troubleShoot = 0; % troubleshooting? y = 1, n = 0

% search volume (note: cylinder separation is 14.5 mm in the phantom, 
% circle diameter is approximately 13 mm)

% determined by the non-upscaled correlation coefficient search to calculate 
% the weigthed sum

% find necessary range to sample
xRange = floor(searchDist/voxelWidth);
yRange = floor(searchDist/voxelHeight);
xRangeWeight = floor(searchDistWeight/voxelWidth);
yRangeWeight = floor(searchDistWeight/voxelHeight);


volXDim = length(imgData(1,:));
volYDim = length(imgData(:,1));

% range of the cylinder (will subtract from the range needed to be tested because
% the cylinder itself occupies some volume within the search area)
cylinderXWidth = ceil(length(circleData(1,:))/2);
cylinderYWidth = ceil(length(circleData(:,1))/2);
cylinderXDim = length(circleData(1,:));
cylinderYDim = length(circleData(:,1));

% % normalize cylinder to signal in the entire region that will be sampled
normX = (startLocation(1) - xRange - cylinderXWidth):1:(startLocation(1) + xRange + cylinderXWidth);
normY = (startLocation(2) - yRange - cylinderYWidth):1:(startLocation(2) + yRange + cylinderYWidth);
% % avoid going beyond bounds of image
if max(normX) > volXDim
    normX = (startLocation(1) - xRange):1:volXDim;
end
if max(normY) > volYDim
    normY = (startLocation(2) - yRange):1:volYDim;
end
if min(normX) < 1
    normX = 1:1:(startLocation(1) + xRange);
end
if min(normY) < 1
    normY = 1:1:(startLocation(2) + yRange);
end

% normalize cylinder to signal in the volume to be searched
% normFactor = max(max(max(volSigData(normY,normX,normZ))));
normFactor = 1; % no normalization to signal in the region
normCircle = circleData.*normFactor;

% the bounds will be determined by the center location, the volume that
% will be searched, and the volume of the cylinder
xMin = (startLocation(1) - xRange - cylinderXWidth);
xMax = (startLocation(1) + xRange - cylinderXWidth);
yMin = (startLocation(2) - yRange - cylinderYWidth);
yMax = (startLocation(2) + yRange - cylinderYWidth);

% make sure to stay within bounds of the image and avoid errors
if xMax > (volXDim - cylinderXDim + 1)
    xMax = volXDim - cylinderXDim + 1; % avoid going out of image volume
end
if xMax < 1
    xMax = 1;
end
if xMin < 1
    xMin = 1; % avoid going out of image volume
end
if yMax > (volYDim - cylinderYDim + 1)
    yMax = volYDim - cylinderYDim + 1; % avoid going out of image volume
end
if yMax < 1
    yMax = 1;
end
if yMin < 1
    yMin = 1; % avoid going out of image volume
end
% make sure xMin not > xMax
if xMin > xMax
    xMin = xMax;
end
if yMin > yMax
    yMin = yMax;
end

% make sure the bounds are integer values
xMin = round(xMin);
xMax = round(xMax);
yMin = round(yMin);
yMax = round(yMax);
% initialize loop
optCorrelation = 0; % optimal correlation
% will return initial guess w/ optCorrelation = 0 if corrcoeff = nan
xCoord = startLocation(1); 
yCoord = startLocation(2);
for X = xMin:1:xMax;
    % x-indices to search
    xSearch = (X:1:(X+cylinderXDim - 1));
    for Y = yMin:1:yMax;
        % y-indices to search
        ySearch = (Y:1:(Y+cylinderYDim - 1));
        % calculate the correlation between the circle and the current
        % image selected in the loop (minus 1 otherwise dimensions
        % will not agree)
        % rotated because of way data is stored
        currentImg = imgData(ySearch,xSearch);
        currentCorrelation = corrcoef(normCircle,currentImg);
        if ((currentCorrelation(2,1)) > optCorrelation)
            % a new optimal correlation, update cylinder location
            optCorrelation = currentCorrelation(2,1);
            % the corrdinates are at the median of the region sampled
            xCoord = median(xSearch);
            yCoord = median(ySearch);
            plotPhantomVol = currentImg; % store the data for plotting
            optXInd = xSearch; % optimal x-indices
            optYInd = ySearch; % optimal y-indices
        end
    end
end
coordCylinderPlot = [yCoord,xCoord]; % The coordinates of the cylinder for plotting 
coordCylinderData = [xCoord,yCoord]; % The coordinates of the cylinder as they appear to user in image

% the weighted sum Method
xMinWeight = (coordCylinderData(1) - xRangeWeight);
if xMinWeight < 1
   xMinWeight = 1; 
end

xMaxWeight = (coordCylinderData(1) + xRangeWeight);
if xMaxWeight > volXDim
    xMaxWeight = volXDim;
end
if xMaxWeight < 1
    xMaxWeight = 1;
end
yMinWeight = (coordCylinderData(2) - yRangeWeight);
if yMinWeight < 1
   yMinWeight = 1; 
end
yMaxWeight = (coordCylinderData(2) + yRangeWeight);
if yMaxWeight > volYDim
    yMaxWeight = volYDim;
end
if yMaxWeight < 1
    yMaxWeight = 1;
end
% make sure min is not greater than max
if xMinWeight > xMaxWeight
    xMinWeight = xMaxWeight;
end
if yMinWeight > yMaxWeight
    yMinWeight = yMaxWeight;
end
xSearchWeight = xMinWeight:1:xMaxWeight;
ySearchWeight = yMinWeight:1:yMaxWeight;
sigImgWeight = imgData(ySearchWeight,xSearchWeight); % reversed data because of way it is stored
weightCoord = weightSumLoc2D(xSearchWeight,ySearchWeight,sigImgWeight);
weightCoordPlot = [weightCoord(2),weightCoord(1)];
% calculate location of center based on upscaled volume data%
%
% NOTE: The input index of the x and y corrdinates for startLocation are
% reversed for this function to findcylinderCenter.m This is due to the way
% that interp3 works and plotting in imshow is indexed.
if upScale == 1;
    [coordCylinderDataUp, optCorrelationUp] = ...
        upscaleSearch2D(imgData,coordCylinderData,voxelHeight,voxelWidth,cylinderBoundary,MRorCT,upscaleFactor);
    % nearst voxel in original image using the upscaled cylinder location
    coordCylinderPlotUp = round([coordCylinderDataUp(2),coordCylinderDataUp(1)]);
else
    coordCylinderDataUp = nan;
    optCorrelationUp = nan;
    coordCylinderPlotUp = nan;
end




%% troubleshooting plots and GUI
if troubleShoot == 1;
    % plot the results
    figure;
    windHigh = normFactor; % normalize windowing
    windLow = 0; % normalize windowing
    subplot(3,1,1)
    imshow(plotPhantomVol(:,:),[windLow windHigh])
    title('plotPhantomVol Image')
    subplot(3,1,2)
    volData = imgData(optYInd,optXInd);
    imshow(volData(:,:),[windLow windHigh])
    title('Volume Data Image')
    subplot(3,1,3)
    imshow(normCircle(:,:),[windLow windHigh])
    title('Circle Image')
    figure;
    imshow(imgData(:,:),[windLow windHigh])
    title('Entire Image','FontSize',24)
end

end