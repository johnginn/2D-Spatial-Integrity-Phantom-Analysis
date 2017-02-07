% This script cycles through the different cylinders in the spatial integrity
% phantom in order to calculate their location.
%
% Input:
% imgData The volume of signal data (dimensions: x,y)
% upscaleImg The upscaled image data (dimensions: x,y)
% circleData The volume of the cylinder (dimensions: x,y)
% voxelHeight (mm/voxel) The height of the pixel
% voxelWidth (mm/voxel) The width of the pixel
% firstCol The index of the first column
% firstRow The index of the first row
% nRows The number of rows
% searchDist (mm) distance in each direction to search for the cylinder
% searchDistWeight (mm) distance in each direction from the location
% cylinderBoundary  y = 1, n = 0 Whether or not a boundary of 0's surrounds the cylinder model
% upScale y = 1, n = 0 search for the cylinder center based off of upscaled data
% MRorCT ('mr' or 'ct') Determines whether the template contrast exhibits a CT or MR scan 
% upscaleFactor The factor the image is upscaled by
% missingLocations The missing cylinder locations [missingRow1,missingCol1;missingRow2,missingCol2;missingRow3,missingCol3];
% 
% Output:
% coordCylinderData [xPos yPos] The position of the cylinder located by program 
% coordcylinderPlot [yPos xPos] The position of the cylinder located by program 
% optCorrelation  The optimal correlation coefficient for each cylinder
% cylGoundTruth [yPos xPos] The starting location of the search area (taken to be
% the ground truth)
% weightCoord [xPos yPos] The position of the cylinder by the weighted sum method
% weightCoordPlot [yPos xPos] The position of the cylinder by the weighted sum method
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
% volSearchRegion Volume of 1's and 0's showing the search region used in
% the analysis
%
% volAllGndTruth Volume of 1's and 0's forming + shape for checking ground
% truth locations
%
% volWeightSum Volume of 1's and 0's forming + shape for checking weighted
% sum calculation regions
%
% John Ginn
% Created: 11/2/16
% Modified: 12/13/16

function [coordCylinderData, coordCylinderPlot, optCorrelation, cylGoundTruth,...
    weightCoord, weightCoordPlot, coordCylDataUp,coordCylPlotUp, optCorrelationUp,...
    volSearchRegion,volAllGndTruth,volWeightSum] = ...
    PhantomScan2D(imgData,upscaleImg,circleData, voxelHeight,voxelWidth,...
    tenthCol,tenthRow,nRows,searchDist,searchDistWeight,cylinderBoundary,MRorCT,upscaleFactor,missingLocations,distBtwnCircles)
% always do the upscaled search
upScale = 1;

% phantom information
nCol = 20; % the number of columns of cylinders
nTotData = 397;
% distBtwnCircles = 14.5; % (mm) The distance between the cylinders in each direction
xIndPerSph = distBtwnCircles/voxelWidth; % number of pixels in x-direction between the cylinders
yIndPerSph = distBtwnCircles/voxelHeight; % number of pixels in y-direction between the cylinders
searchDistStore = searchDist; % store the search distance so that we can modify it for the spheres in the bottom right hand corner


% determine the progress of finding the cylinders
percentFound = 10:10:100;
numcylindersPercentFound = floor(nTotData./10.*(1:10)); % array of number of cylinders corresponding to percent found
countPercentFound = 1;

% initialize data for efficiency
coordCylinderData = cell(1,nTotData);
coordCylinderPlot = cell(1,nTotData);
optCorrelation = cell(1,nTotData);
cylGoundTruth = cell(1,nTotData);
coordCylDataUp = cell(1,nTotData);
coordCylPlotUp = cell(1,nTotData);
optCorrelationUp = cell(1,nTotData);
weightCoord  = cell(1,nTotData);
weightCoordPlot  = cell(1,nTotData);
searchIndX = cell(1,nTotData);
searchIndY = cell(1,nTotData);
% array of the weighted sum distances used to calculate the locations
currWeightSumRange = searchDistWeight.*ones(1,nTotData);
% cycle through slices
saveCount = 1; % used to save the data
% step through the rows
for rowStep = 1:nRows
    % the count should start at 0 so that the position does not move at the
    % first iteration
    rowPos = rowStep - 10;
    currentRow = tenthRow + round(yIndPerSph*rowPos);
    gndTruthRow = tenthRow + yIndPerSph*rowPos;
    for colStep = 1:nCol
        % the count should start at 0 so that the position does not move
        % at the first iteration
        colPos = colStep - 10; % remove one since the first step
        % was used to calculate the center location
        % the current column location
        currentCol = tenthCol + round(xIndPerSph*colPos);
        gndTruthCol = tenthCol + xIndPerSph*colPos;
        % Avoid the three "missing" cylinders
        % real time scan
        if ((rowStep == missingLocations(1,1))&&(colStep == missingLocations(1,2)))||...
                ((rowStep == missingLocations(2,1))&&(colStep == missingLocations(2,2)))||...
                ((rowStep == missingLocations(3,1))&&(colStep == missingLocations(3,2)))
            % missing cylinder, do nothing
        else
            startLocation = [currentCol currentRow];
            if (rowStep > nRows - 3)&&(colStep == nCol)
                searchDist = 0; % 
            else
                searchDist = searchDistStore;
            end
            % Find the cylinder
            [coordCylinderData{saveCount},coordCylinderPlot{saveCount}, optCorrelation{saveCount},...
                weightCoord{saveCount}, weightCoordPlot{saveCount},...
                coordCylDataUp{saveCount}, coordCylPlotUp{saveCount}, optCorrelationUp{saveCount},...
                searchIndX{saveCount},searchIndY{saveCount}] = ...
                findCircleCenter2D(imgData, circleData,startLocation,voxelHeight,voxelWidth,...
                searchDist,currWeightSumRange(saveCount),upScale,cylinderBoundary,MRorCT,upscaleFactor);
            % print progress of cycling through the cylinders
            if saveCount == numcylindersPercentFound(countPercentFound)
                disp(strcat([num2str(percentFound(countPercentFound)),...
                    '% of cylinder locations determined']))
                countPercentFound = countPercentFound + 1;
            end
            % ground truth locations (note rotation for plotting) [y x z]
            cylGoundTruth{saveCount} = [gndTruthRow gndTruthCol];
            saveCount = saveCount + 1; % used to save the data
        end
    end
end


% initialize volume of data for search region using correlation coeff
volDimX = length(imgData(:,1));
volDimY = length(imgData(1,:));
volSearchRegion = zeros(volDimX,volDimY);
volAllGndTruth = volSearchRegion;
volWeightSum = volSearchRegion;
% step through each of the cylinders and construct a border to show search
% region
for step = 1:(saveCount - 1)
    % make volume of + shapes for plotting
    gndTruthCurrent = round(cylGoundTruth{step});
    gndTruthX = (gndTruthCurrent(1)-1):1:(gndTruthCurrent(1)+1);
    gndTruthY = (gndTruthCurrent(2)-1):1:(gndTruthCurrent(2)+1);
    % no need to rotate x,y because that was already done
    volAllGndTruth(gndTruthX,gndTruthCurrent(2)) = 1;
    volAllGndTruth(gndTruthCurrent(1),gndTruthY) = 1;
    % search region by correlation coefficent
    xRegion = searchIndX{step};
    xMin = min(xRegion);
    xMax = max(xRegion);
    yRegion = searchIndY{step};
    yMin = min(yRegion);
    yMax = max(yRegion);
    for stepX = xRegion
        for stepY = yRegion
            % if at one of the boundaries of the search region store a
            % 1 for plotting
            if (((stepX == xMin) || (stepX == xMax))||...
                    ((stepY == yMin) || (stepY == yMax)))
                % ensure integer values
                stepX = round(stepX);
                stepY = round(stepY);
                % reverse for plotting
                volSearchRegion(stepY,stepX) = 1;
            end
        end
    end
    
    % the location found by the correlation coefficent for the current
    % cylinder
    cylinderCurrent = round(coordCylinderPlot{step});
    % weighted sum region calculation
    xWeightDist = round(currWeightSumRange(step)/voxelWidth);
    yWeightDist = round(currWeightSumRange(step)/voxelHeight);
    % x-bounds
    if (cylinderCurrent(1)-xWeightDist) < 1;
        xMinWeight = 1;
    else
        xMinWeight = (cylinderCurrent(1)-xWeightDist);
    end
    if (cylinderCurrent(1)+xWeightDist) > volDimX;
        xMaxWeight = xMax;
    else
        xMaxWeight = (cylinderCurrent(1)+xWeightDist); % already ensured this does not go beyond boundary
    end
    % y-bounds
    if (cylinderCurrent(2)-yWeightDist) < 1;
        yMinWeight = 1;
    else
        yMinWeight = (cylinderCurrent(2)-yWeightDist);
    end
    if (cylinderCurrent(2)+yWeightDist) > volDimY;
        yMaxWeight = yMax;
    else
        yMaxWeight = (cylinderCurrent(2)+yWeightDist); % already ensured this does not go beyond boundary
    end
    
    % make sure min is not > max
    if xMinWeight > xMaxWeight
        xMinWeight = xMaxWeight;
    end
    if yMinWeight > yMaxWeight
        yMinWeight = yMaxWeight;
    end
    xWeight = xMinWeight:1:xMaxWeight;
    yWeight = yMinWeight:1:yMaxWeight;
    for stepX = xWeight
        for stepY = yWeight
            % if at one of the boundaries of the search region store a
            % 1 for plotting
            if (((stepX == xMinWeight) || (stepX == xMaxWeight))||...
                    ((stepY == yMinWeight) || (stepY == yMaxWeight)))
                % ensure integer values
                stepX = round(stepX);
                stepY = round(stepY);
                volWeightSum(stepX,stepY) = 1;
            end
        end
    end
end


end