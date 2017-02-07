% Function to compare this analysis to the analysis code provided by
% ViewRay using the center transverse slice data.
%
% Input:
% slicecylinderX The x-location of the cylinders in the slice extracted
% slicecylinderY The y-location of the cylinders in the slice extracted
% sliceGroundX The x-location of the ground truth in the slices
% sliceGroundY The y-location of the ground truth in the slices
% voxelHeight (mm) The height of the voxels
% voxelWidth (mm) The width of the voxels
% centSliceImg The image data for the phantom
% radiusThreshold1 The deviation tolerance first threshold to calculate the percent passing 
% radiusThreshold2 The deviation tolerance second threshold to calculate the percent passing 
% radiusThreshold3 The deviation tolerance third threshold to calculate the percent passing 
% centerCol The index of the center column of the phantom
% centerRow The index of the center row of the phantom
% radiusSearch1 (mm) The radius from isocenter that defines the subvolume for the first threshold
% radiusSearch2 (mm) The radius from isocenter that defines the subvolume for the second threshold
% radiusSearch3 (mm) The radius from isocenter that defines the subvolume for the third threshold
% scaleFactor The factor that the image is upscaled by
% sliceShift The distance away from isocenter the phantom was shifted during scanning
% xCircle The x-location of the cylinder found before any rotation corrections
% yCircle The y-location of the cylinder found before any rotation corrections
% xCorrGndTruth The x-ground truth location after rotation correction
% yCorrGndTruth The y-ground truth location after rotation correction
% plotData Whether or not to plot the data
% plotOutsideRegion Whether or not to plot the cylinders passing or failing radiusThershol3 outside radius 3
% MRorCT Whether the scan is an 'mr' scan or a 'ct' scan
% 
% Output:
% compareDeviationXY The deviation between each cylinder location
% and ground truth in the X,Y plane
% avgDevCompViewRayXY The average deviation between the cylinder location
% and ground truth in the X,Y plane
% percentPassXY1 The percentage points passing only considering deviation in two dimensions (1st threshold)
% percentPassXY2 The percentage points passing only considering deviation in two dimensions (2nd threshold)
% percentPassXY3 The percentage points passing only considering deviation in two dimensions (3rd threshold)
% avgDevRadius1XY Average 2D deviation for cylinders within radius 1
% avgDevRadius2XY Average 2D deviation for cylinders within radius 2
% avgDevRadius3XY Average 2D deviation for cylinders within radius 3
% maxDevRadius1XY Maximum 2D deviation for cylinders within radius 1
% maxDevRadius2XY Maximum 2D deviation for cylinders within radius 2
% maxDevRadius3XY Maximum 2D deviation for cylinders within radius 3
%
% John Ginn
% Created: 11/3/16
% Modified: 11/4/16
function [compareDeviationXY,avgDevAllXY,percentPassXY1,percentPassXY2,percentPassXY3,...
    avgDevRadius1XY,avgDevRadius2XY,avgDevRadius3XY,...
     maxDevRadius1XY,maxDevRadius2XY,maxDevRadius3XY] = ...
    compViewRay2D(slicecylinderX,slicecylinderY,sliceGroundX,sliceGroundY,...
    voxelHeight,voxelWidth,centSliceImg,radiusThreshold1,radiusThreshold2,radiusThreshold3,...
    centerCol,centerRow,radiusSearch1,radiusSearch2,radiusSearch3,scaleFactor,sliceShift,...
    xCircle,yCircle,xCorrGndTruth,yCorrGndTruth,plotData,plotOutsideRegion,MRorCT)
%% scale factor for upscaling the image
% calcualte the deviation for the center transverse slice to compare to the
% ViewRay analysis software
avgDevAllXY = 0;
compareDeviationXY = zeros(1,length(slicecylinderX));
radius = zeros(1,length(slicecylinderX));
if (~strcmp(MRorCT,'ct'))
    cylinderFitPtsPass = zeros(scaleFactor*length(centSliceImg(:,1)),scaleFactor*length(centSliceImg(1,:)));
else
    cylinderFitPtsPass = zeros(length(centSliceImg(:,1)),length(centSliceImg(1,:)));
end
cylinderFitPtsFail = cylinderFitPtsPass;
cylinderPtsOutsideRegionsPass = cylinderFitPtsPass;
cylinderPtsOutsideRegionsFail = cylinderFitPtsPass;
grndTruthPts = cylinderFitPtsPass;
countPassXY1 = 0; % 1st threshold
countPassXY2 = 0; % 2nd threshold
countPassXY3 = 0; % 3nd threshold
countcylindersInRadius1 = 0; % the number of cylinders within the 1st specified radius 
countcylindersInRadius2 = 0; % the number of cylinders within the 2nd specified radius 
countcylindersInRadius3 = 0; % the number of cylinders within the 3rd specified radius 
% relative to isocenter
% specifiy not a number in case sliceShift > any of the radius of the search regions
devSmallXY1 = nan;
devSmallXY2 = nan;
devSmallXY3 = nan;
percentPassXY1 = nan;
percentPassXY2 = nan;
percentPassXY3 = nan;
for step = 1:length(slicecylinderX)    
    % current distance from isocenter
    radius(step) = sqrt((voxelWidth*(slicecylinderX(step) - centerCol)).^2 + ...
        (voxelHeight*(slicecylinderY(step) - centerRow)).^2 +...
        sliceShift.^2);
    
    % calculate x,y deviation
    compareDeviationXY(step) =  sqrt((voxelWidth*(slicecylinderX(step) - sliceGroundX(step))).^2 + ...
        (voxelHeight*(slicecylinderY(step) - sliceGroundY(step))).^2);
    avgDevAllXY = avgDevAllXY + compareDeviationXY(step);

    % check to see if point is within specified distance of isocenter
    if radius(step) <= radiusSearch1
        countcylindersInRadius1 = countcylindersInRadius1 + 1;
        radiusSmall1(countcylindersInRadius1) = radius(step);
        devSmallXY1(countcylindersInRadius1) = compareDeviationXY(step);
        % check XY deviation to see if it passes the threshold
        if compareDeviationXY(step) <= radiusThreshold1
            countPassXY1 = countPassXY1 + 1;
        end
    end
    if radius(step) <= radiusSearch2
        countcylindersInRadius2 = countcylindersInRadius2 + 1;
        radiusSmall2(countcylindersInRadius2) = radius(step);
        devSmallXY2(countcylindersInRadius2) = compareDeviationXY(step);
        % check XY deviation to see if it passes the threshold
        if compareDeviationXY(step) <= radiusThreshold2
            countPassXY2 = countPassXY2 + 1;
        end
    end
    if radius(step) <= radiusSearch3
        countcylindersInRadius3 = countcylindersInRadius3 + 1;
        radiusSmall3(countcylindersInRadius3) = radius(step);
        devSmallXY3(countcylindersInRadius3) = compareDeviationXY(step);
        % check XY deviation to see if it passes the threshold
        if compareDeviationXY(step) <= radiusThreshold3
            countPassXY3 = countPassXY3 + 1;
        end
    end
    % round locations for plotting, reverse points because of imshow
    % indexing
%     data = round([slicecylinderY(step) slicecylinderX(step)].*scaleFactor);
    % use the non-rotated found location
    if (~strcmp(MRorCT,'ct'))
        data = round([yCircle(step) xCircle(step)].*scaleFactor);
        xMarkerLoc = (data(1)-scaleFactor - 1):1:(data(1)+scaleFactor + 1); % make + symbol on image
        yMarkerLoc = (data(2)-scaleFactor - 1):1:(data(2)+scaleFactor + 1);% make + symbol on image
    else
        data = round([yCircle(step) xCircle(step)]);
        xMarkerLoc = (data(1)-2):1:(data(1)+2); % make + symbol on image
        yMarkerLoc = (data(2)-2):1:(data(2)+2);% make + symbol on image
    end
    
    % choose dataset for determining plot pass/fail
    plotDev = compareDeviationXY(step);
    plotRadius = radius(step);
    
    % first determine the location of the cylinder
    if (plotRadius < radiusSearch1)
        if plotDev <= radiusThreshold1
            % the marker for the cylinders that passed
            cylinderFitPtsPass(xMarkerLoc,yMarkerLoc) = makeShape(length(xMarkerLoc),'square','thick'); % make a square
        else % the marker for the cylinders that failed
            cylinderFitPtsFail(xMarkerLoc,yMarkerLoc) = makeShape(length(xMarkerLoc),'x','thick'); % make a x
        end
    elseif (plotRadius >= radiusSearch1)&&(plotRadius < radiusSearch2)
        if plotDev <= radiusThreshold2
            % the marker for the cylinders that passed
            cylinderFitPtsPass(xMarkerLoc,yMarkerLoc) = makeShape(length(xMarkerLoc),'square','thick'); % make a square
        else % the marker for the cylinders that failed
            cylinderFitPtsFail(xMarkerLoc,yMarkerLoc) = makeShape(length(xMarkerLoc),'x','thick'); % make a x
        end
    elseif (plotRadius >= radiusSearch2)&&(plotRadius < radiusSearch3)
        if plotDev <= radiusThreshold3
            % the marker for the cylinders that passed
            cylinderFitPtsPass(xMarkerLoc,yMarkerLoc) = makeShape(length(xMarkerLoc),'square','thick'); % make a square
        else % the marker for the cylinders that failed
            cylinderFitPtsFail(xMarkerLoc,yMarkerLoc) = makeShape(length(xMarkerLoc),'x','thick'); % make a x
        end
    end
    % plot regions outside the three typical regions for debugging
    if plotOutsideRegion == 1;
        if (plotRadius >= radiusSearch3)
            if plotDev <= radiusThreshold3
                % the marker for the cylinders that passed
                cylinderFitPtsPass(xMarkerLoc,yMarkerLoc) = makeShape(length(xMarkerLoc),'square','thick'); % make a square
            else % the marker for the cylinders that failed
                cylinderFitPtsFail(xMarkerLoc,yMarkerLoc) = makeShape(length(xMarkerLoc),'x','thick'); % make a x
            end
        end
    end

    % Marker for fusing on the image (add on data in each slice)
    % marker{sliceNumber} = [current data; new data]
    % store location of ground truth for cylinders
%     data = round([sliceGroundY(step) sliceGroundX(step)].*scaleFactor);
    % use the rotated ground truth data
    if (~strcmp(MRorCT,'ct'))
        data = round([yCorrGndTruth(step) xCorrGndTruth(step)].*scaleFactor);
        xMarkerLoc = (data(1)-scaleFactor):1:(data(1)+scaleFactor); % make + symbol on image
        yMarkerLoc = (data(2)-scaleFactor):1:(data(2)+scaleFactor);% make + symbol on image
    else
        data = round([yCorrGndTruth(step) xCorrGndTruth(step)]);
        xMarkerLoc = (data(1)-1):1:(data(1)+1); % make + symbol on image
        yMarkerLoc = (data(2)-1):1:(data(2)+1);% make + symbol on image
    end
    if plotOutsideRegion ~= 1;
        if (plotRadius <= radiusSearch3)
            grndTruthPts(xMarkerLoc,data(2)) = 1; % - portion of + marker for ground truth of cylinders
            grndTruthPts(data(1),yMarkerLoc) = 1; % - portion of + marker for ground truth of cylinders
            grndTruthPts(xMarkerLoc,data(2)-1) = 1; % - portion of + marker for ground truth of cylinders
            grndTruthPts(data(1)-1,yMarkerLoc) = 1; % - portion of + marker for ground truth of cylinders
            grndTruthPts(xMarkerLoc,data(2)+1) = 1; % - portion of + marker for ground truth of cylinders
            grndTruthPts(data(1)+1,yMarkerLoc) = 1; % - portion of + marker for ground truth of cylinders
        end
    else
        grndTruthPts(xMarkerLoc,data(2)) = 1; % - portion of + marker for ground truth of cylinders
        grndTruthPts(data(1),yMarkerLoc) = 1; % - portion of + marker for ground truth of cylinders
        grndTruthPts(xMarkerLoc,data(2)-1) = 1; % - portion of + marker for ground truth of cylinders
        grndTruthPts(data(1)-1,yMarkerLoc) = 1; % - portion of + marker for ground truth of cylinders
        grndTruthPts(xMarkerLoc,data(2)+1) = 1; % - portion of + marker for ground truth of cylinders
        grndTruthPts(data(1)+1,yMarkerLoc) = 1; % - portion of + marker for ground truth of cylinders
    end

end
avgDevRadius1XY = sum(devSmallXY1)/length(devSmallXY1);
avgDevRadius2XY = sum(devSmallXY2)/length(devSmallXY2);
avgDevRadius3XY = sum(devSmallXY3)/length(devSmallXY3);
maxDevRadius1XY = max(devSmallXY1);
maxDevRadius2XY = max(devSmallXY2);
maxDevRadius3XY = max(devSmallXY3);
avgDevAllXY = avgDevAllXY/length(slicecylinderX);
percentPassXY1 = 100*countPassXY1/countcylindersInRadius1;
percentPassXY2 = 100*countPassXY2/countcylindersInRadius2;
percentPassXY3 = 100*countPassXY3/countcylindersInRadius3;

if (~strcmp(MRorCT,'ct'))
    % upscale the image for plotting
    % original information
    origX = 1:length(centSliceImg(1,:));
    origY = 1:length(centSliceImg(:,1));
    [origMeshX,origMeshY] = meshgrid(origX,origY);
    % new data
    volUpscaleX = 1/scaleFactor:1/scaleFactor:length(centSliceImg(1,:));
    volUpscaleY = 1/scaleFactor:1/scaleFactor:length(centSliceImg(:,1));
    [upMeshX, upMeshY] = meshgrid(volUpscaleX,volUpscaleY);
    % the upscaled data
    SigDataUpscale = interp2(origMeshX,origMeshY,centSliceImg,...
        upMeshX,upMeshY,'cubic');
    
    upVoxWidth = voxelWidth/scaleFactor;
    upVoxHeight = voxelHeight/scaleFactor;
    % check if dimensions are even or odd. If odd, just divide by 2 (center
    % position is an integer value), otherwise add 0.5
    if mod(length(SigDataUpscale(:,1))/2,2) ~= 0
        upCenterRow = length(SigDataUpscale(:,1))/2;
    else
        upCenterRow = length(SigDataUpscale(:,1))/2 + 0.5;
    end
    if mod(length(SigDataUpscale(1,:))/2,2) ~= 0
        upCenterCol = length(SigDataUpscale(1,:))/2;
    else
        upCenterCol = length(SigDataUpscale(1,:))/2 + 0.5;
    end
else
    SigDataUpscale = centSliceImg;
    upCenterRow = centerRow;
    upCenterCol = centerCol;
    upVoxWidth = voxelWidth;
    upVoxHeight = voxelHeight;
end

if plotData == 1
    figure;
    % fit marker plotting
    imshow(SigDataUpscale,[])
    % search region boundary plotting
    [circleImg1] = imageCircle(SigDataUpscale,radiusSearch1,upVoxHeight,upVoxWidth,...
        upCenterRow,upCenterCol,sliceShift);
    [circleImg2] = imageCircle(SigDataUpscale,radiusSearch2,upVoxHeight,upVoxWidth,...
        upCenterRow,upCenterCol,sliceShift);
    [circleImg3] = imageCircle(SigDataUpscale,radiusSearch3,upVoxHeight,upVoxWidth,...
        upCenterRow,upCenterCol,sliceShift);
    for step = 1:3
        if step == 1
            currentImg = circleImg1;
        elseif step == 2
            currentImg = circleImg2;
        else
            currentImg = circleImg3;
        end
        rgb= [0,183,229]/255;
        blue = cat(3, rgb(1).*ones(size(currentImg)),rgb(2).*ones(size(currentImg)),rgb(3).*ones(size(currentImg)));
        hold on
        hBlue = imshow(blue);
        hold off
        set(hBlue, 'AlphaData', currentImg) % make color sheet only show markers
    end
    
    % cat(3,r,g,b)
    green = cat(3, zeros(size(cylinderFitPtsPass)),ones(size(cylinderFitPtsPass)), zeros(size(cylinderFitPtsPass)));
    % bright = 0.9;
    % white version
    % green = bright.*cat(3, ones(size(cylinderFitPtsPass)),ones(size(cylinderFitPtsPass)), ones(size(cylinderFitPtsPass)));
    hold on
    hGreen = imshow(green);
    hold off
    set(hGreen, 'AlphaData', cylinderFitPtsPass) % make color sheet only show markers
    % the cylinders that failed
    red = cat(3, ones(size(cylinderFitPtsFail)),zeros(size(cylinderFitPtsFail)), zeros(size(cylinderFitPtsFail)));
    % white version
    % red = bright.*cat(3, ones(size(cylinderFitPtsFail)),ones(size(cylinderFitPtsFail)), ones(size(cylinderFitPtsFail)));
    hold on
    hRed = imshow(red);
    hold off
    set(hRed, 'AlphaData', cylinderFitPtsFail) % make color sheet only show markers
    
    % ground truth marker plotting
    % cat(3,r,g,b)
    white = cat(3, ones(size(grndTruthPts)),ones(size(grndTruthPts)),ones(size(grndTruthPts)));
    hold on
    hWhite = imshow(white);
    hold off
    set(hWhite, 'AlphaData', grndTruthPts) % make color sheet only show markers
    % title('Spatial Integrity Phantom Center Slice','FontSize',20)
    xlabel('y-position (mm)','FontSize',20)
    ylabel('z-position (mm)','FontSize',20)
    title(strcat(['Cylinder Deviation Test, ',num2str(sliceShift),' (mm) shift']),'FontSize',20)
    
    
    % custom axes to show distances (the locations of the axes)
    numOfTicks = 10;
    xAxisSpacing = floor(length(SigDataUpscale(1,:))/numOfTicks);
    yAxisSpacing = floor(length(SigDataUpscale(:,1))/numOfTicks);
    
    % xAxisLocation = linspace(1,length(SigDataUpscale(1,:)),numOfTicks);
    % yAxisLocation = linspace(1,length(SigDataUpscale(:,1)),numOfTicks);
    xAxisLocation = 1:xAxisSpacing:(xAxisSpacing*numOfTicks+1);
    yAxisLocation = 1:yAxisSpacing:(yAxisSpacing*numOfTicks+1);
    % the values on the axes
    if (~strcmp(MRorCT,'ct'))
        xAxisValue = voxelWidth/scaleFactor.*(xAxisLocation - round(median(xAxisLocation)));
        yAxisValue = voxelHeight/scaleFactor.*(yAxisLocation - round(median(yAxisLocation)));
    else
        xAxisValue = voxelWidth.*(xAxisLocation - round(median(xAxisLocation)));
        yAxisValue = voxelHeight.*(yAxisLocation - round(median(yAxisLocation)));
    end
    yAxisLocation(length(yAxisLocation)) = length(SigDataUpscale(:,1));
    xAxisLocation(length(xAxisLocation)) = length(SigDataUpscale(1,:));
    axis on
    axisHandle = gca;
    set(gca,'XTickMode','manual')
    set(gca,'YTickMode','manual')
    set(gca,'XTick',xAxisLocation)
    set(gca,'YTick',yAxisLocation)
    set(gca,'XTickLabel',xAxisValue)
    set(gca,'YTickLabel',yAxisValue)
    set(gca,'FontSize',16)
    
end


end