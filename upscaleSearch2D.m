% This function takes the locations of the cylinders that were found using
% the native resolution images, upscales the volume data and extracts a
% more precise location of the cylinders.
%
% Input:
% volSigData The volume of signal data (dimensions: x,y)
% centerSearch [yInd,xInd] The center of the search volume
% voxelHeight (mm/voxel) The height of the pixel
% voxelWidth (mm/voxel) The width of the pixel
% cylinderBoundary  y = 1, n = 0 Whether or not a boundary of 0's surrounds the cylinder model
% MRorCT ('mr' or 'ct') Determines whether the template contrast exhibits a CT or MR scan 
% upscaleFactor The factor the image is upscaled by
%
% Output:
% coordCylinderDataUp [xInd,yInd] The coordinates of the cylinder from upscaled
% analysis as the appear to user in imshow
%
% optCorrelationUp The optimal correlation coefficient from upscaled
% analysis
%
% NOTE: The input index of the x and y corrdinates for centerSearch are
% reversed for this function to findcylinderCenter.m This is due to the way
% that interp3 works and plotting in imshow is indexed.
%
% John Ginn
% Created: 11/2/16
% Modified: 11/4/16

function [coordCylinderDataUp, optCorrelationUp] =...
    upscaleSearch2D(imgSigData,centerOfCirc,voxelHeight,voxelWidth,cylinderBoundary,MRorCT,upscaleFactor)
troubleShoot = 0; % troubleshooting? y = 1, n = 0

circleDiameter = 13; % (mm)
searchDist = 3; % (mm) distance in each direction to search for the cylinder

upscaleVoxX = voxelWidth/upscaleFactor; % (mm) the desired upscale voxel x-dimension
upscaleVoxY = voxelHeight/upscaleFactor; % (mm) the desired upscale voxel y-dimension
scaleFactorX = upscaleFactor;
scaleFactorY = upscaleFactor;


% calculate dimensions of upscaled data (calculating entire volume would
% take too long and require more memory)
volXDimOrig = length(imgSigData(1,:)); % before scaling
volYDimOrig = length(imgSigData(:,1)); % before scaling
volXDim = scaleFactorX.*volXDimOrig; % after scaling
volYDim = scaleFactorY.*volYDimOrig; % after scaling

% find necessary range to sample
xRangeOrig = floor(searchDist/voxelWidth); % before scaling
yRangeOrig = floor(searchDist/voxelHeight); % before scaling
xRange = floor(searchDist/upscaleVoxX); % after scaling
yRange = floor(searchDist/upscaleVoxY); % after scaling

% the upscaled cylinder volume makecylinder(height,width,length,diameter)
[circleImgUpscale] = ...
    makeCircle(upscaleVoxY, upscaleVoxX, circleDiameter,cylinderBoundary,MRorCT);
[circleImgOrig] = ...
    makeCircle(voxelHeight, voxelWidth, circleDiameter,cylinderBoundary,MRorCT);

% range of the cylinder (will subtract from the range needed to be tested because
% the cylinder itself occupies some volume within the search area)
% upscale
cylinderXWidth = ceil(length(circleImgUpscale(1,:))/2);
cylinderYWidth = ceil(length(circleImgUpscale(:,1))/2);
cylinderXDim = length(circleImgUpscale(1,:));
cylinderYDim = length(circleImgUpscale(:,1));
% original
cylinderXWidthOrig = ceil(length(circleImgOrig(1,:))/2);
cylinderYWidthOrig = ceil(length(circleImgOrig(:,1))/2);
cylinderXDimOrig = length(circleImgOrig(1,:));
cylinderYDimOrig = length(circleImgOrig(:,1));

% the location of the upscaled center search region
upscCenterSearch(1) = centerOfCirc(1).*scaleFactorX;
upscCenterSearch(2) = centerOfCirc(2).*scaleFactorY;

% normalize cylinder to signal in the entire region that will be sampled
normX = (centerOfCirc(1) - xRangeOrig - cylinderXWidthOrig):1:...
    (centerOfCirc(1) + xRangeOrig + cylinderXWidthOrig);
normY = (centerOfCirc(2) - yRangeOrig - cylinderYWidthOrig):1:...
    (centerOfCirc(2) + yRangeOrig + cylinderYWidthOrig);
% avoid going beyond bounds of image
if max(normX) > volXDimOrig
    normX = (centerOfCirc(1) - xRangeOrig):1:volXDimOrig;
end
if max(normY) > volYDimOrig
    normY = (centerOfCirc(2) - yRangeOrig):1:volYDimOrig;
end
if min(normX) < 1
    normX = 1:1:(centerOfCirc(1) + xRangeOrig);
end
if min(normY) < 1
    normY = 1:1:(centerOfCirc(2) + yRangeOrig);
end
% normalize cylinder to signal in the volume to be searched
normFactor = max(max(imgSigData(normY,normX)));
normcylinder = circleImgUpscale.*normFactor;

% the bounds will be determined by the center location, the volume that
% will be searched, and the volume of the cylinder
xMin = (upscCenterSearch(1) - xRange - cylinderXWidth);
xMax = (upscCenterSearch(1) + xRange - cylinderXWidth);
yMin = (upscCenterSearch(2) - yRange - cylinderYWidth);
yMax = (upscCenterSearch(2) + yRange - cylinderYWidth);

% make sure to stay within bounds of the image and avoid errors
if xMax > (volXDim - cylinderXDim + 1)
    xMax = volXDim - cylinderXDim + 1; % avoid going out of image volume
end
if xMin < 1
    xMin = 1; % avoid going out of image volume
elseif xMin > xMax
    xMin = xMax;
end
if yMax > (volYDim - cylinderYDim + 1)
    yMax = volYDim - cylinderYDim + 1; % avoid going out of image volume
end
if yMin < 1
    yMin = 1; % avoid going out of image volume
elseif yMin > yMax
    yMin = yMax;
end


% calculate the upscaled volume data that will be searched
%
% these bounds are the bounds of the whole subvolume, need to obtain
% addtional region to account for width of cylinder
xMinOrig = (centerOfCirc(1) - xRangeOrig - cylinderXWidthOrig); 
xMaxOrig = (centerOfCirc(1) + xRangeOrig + cylinderXWidthOrig);
yMinOrig = (centerOfCirc(2) - yRangeOrig - cylinderYWidthOrig);
yMaxOrig = (centerOfCirc(2) + yRangeOrig + cylinderYWidthOrig);
% avoid errors going outside the bounds
if xMaxOrig > (volXDimOrig)
    xMaxOrig = volXDimOrig; % avoid going out of image volume
end
if xMinOrig < 1
    xMinOrig = 1; % avoid going out of image volume
elseif xMinOrig > xMaxOrig
    xMinOrig = xMaxOrig;
end
if yMaxOrig > (volYDimOrig)
    yMaxOrig = volYDimOrig; % avoid going out of image volume
end
if yMinOrig < 1
    yMinOrig = 1; % avoid going out of image volume
elseif yMinOrig > yMaxOrig
    yMinOrig = yMaxOrig;
end


%% Upscale the data
origX = xMinOrig:1:xMaxOrig;
origY = yMinOrig:1:yMaxOrig;
[origMeshX,origMeshY] = meshgrid(origX,origY);

% % select out data to be upscaled (store y,x,z like in image data)
imgSigOrig = zeros(length(origY),length(origX));
countX = 1;
countY = 1;
for stepX = origX
    countY = 1; % reset the count
    for stepY = origY       
        countZ = 1; % reset the count
            % the selected data
            imgSigOrig(countY,countX) = imgSigData(stepY,stepX);
        countY = countY + 1;
    end
    countX = countX + 1;
end

% note: do not start at 1/scaleFactor because original data starts at 1 
% starting less than 1 would sample outside of bounds of the volume selected
%
% these are the locations of the interpolation, NOT indexes thus the range
% will be larger than expected by a factor of the width of the cylinder
volUpscaleX = xMinOrig:1/scaleFactorX:xMaxOrig;
volUpscaleY = yMinOrig:1/scaleFactorY:yMaxOrig;
[upMeshX, upMeshY] = meshgrid(volUpscaleX,volUpscaleY);
% the upscaled data
imgSigDataUpscale = interp2(origMeshX,origMeshY,imgSigOrig,...
    upMeshX,upMeshY,'cubic');


% cycle through the volume and calculate correlation coefficients -->
% start beginning correlation as max correlation --> cycle through and
% replace correlation if corrcoef yeilds a higher value

% initialize loop
optCorrelationUp = 0; % optimal correlation
coordCylinderDataUp = centerOfCirc; % optimal coordinates
currentImg = zeros(cylinderXDim,cylinderYDim);
% will return initial guess w/ optCorrelation = 0 if corrcoeff = nan
xCoord = centerOfCirc(1); 
yCoord = centerOfCirc(2);
% plus one on bounds because otherwise last volume location will be skipped
for X = 1:1:(length(imgSigDataUpscale(1,:)) - cylinderXDim + 1);
    % x-indices to search
    xSearch = (X:1:(X+cylinderXDim - 1));
    for Y = 1:1:(length(imgSigDataUpscale(:,1)) - cylinderYDim + 1);
        % y-indices to search
        ySearch = (Y:1:(Y+cylinderYDim - 1));
            % calculate the correlation between the cylinder and the current
            % volume selected in the loop (minus 1 otherwise dimensions
            % will not agree)
            currentImg = imgSigDataUpscale(ySearch,xSearch);
            currentCorrelation = corrcoef(normcylinder,currentImg);
            if ((currentCorrelation(2,1)) > optCorrelationUp)
                % a new optimal correlation, update cylinder location
                optCorrelationUp = currentCorrelation(2,1);
                % the corrdinates are at the median of the region sampled
                xCoordIndex = median(xSearch);
                yCoordIndex = median(ySearch); 
                % because the upscaled sub-volume indexing is different from
                % original volume's index, grab optimal location from meshgrid
                % Note that coordinates are reversed because of meshgrid
                xCoord = upMeshX(1,xCoordIndex);
                yCoord = upMeshY(yCoordIndex,1);
            end
    end
end
% reversed compared to findcylinderCenter
coordCylinderDataUp = [yCoord,xCoord]; % The coordinates of the cylinder as they appear to user in image



%% troubleshooting plots and GUI
if troubleShoot == 1;
    % check the upscale algorithm
    guiData{1} = imgSigOrig;
    guiData{2} = imgSigDataUpscale;
    UpscaleGUI(guiData)
    % plot the results
end
end