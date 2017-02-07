% This function plots the outliers that the user is going to remove in red
% as well as shows the volumetric dataset so that the locations of theses
% cylinders can be confirmed in the phantom volume.
%
% Input:
% xCircleAll (index) The found x-location of the cylinders before corrections are applied in each scan (gnd corrected data)
% yCircleAll (index) The found y-location of the cylinders before corrections are applied in each scan (gnd corrected data)
% finalCircleDataX (index) The found x-location of the cylinders after corrections are applied in each scan (cylinder corrected data)
% finalCircleDataX (index) The found y-location of the cylinders after corrections are applied in each scan (cylinder corrected data)
% deviation2DScans The 2D deviation for each of the cylinders in each scan
% xCorrGndTruthAll (index) The ground-truth x-position after corrections for each scan (gnd corrected data)
% yCorrGndTruthAll (index) The ground-truth y-position after corrections for each scan (gnd corrected data)
% xGndTruthFinalAll (index) The ground-truth x-position before corrections for each scan (cylinder corrected data)
% yGndTruthFinalAll (index) The ground-truth y-position before corrections for each scan (cylinder corrected data)
% radiusAllScans (mm) The distance of the cylinders from isocenter for each scan
% imageDataAll The image of the phantom for each scan
% radiusRemove [min max] (mm) The radii specifying the region of cylinders to be removed
% rangeDevRemove [min max] (mm) The deviations specifying the region of cylinders to be removed
% sliceShiftExcel (mm) The series of slice shifts for each scan
% voxelWidthExcel (mm) The voxel width for each of the scans
% voxelHeightExcel (mm) The voxel height for each of the scans
% radiusSearch1Excel (mm) The radius from isocenter defining the first analysis region for all the scans
% radiusSearch2Excel (mm) The radius from isocenter defining the second analysis region for all the scans
% radiusSearch3Excel (mm) The radius from isocenter defining the third analysis region for all the scans
% scaleFactor The factor the images will be upscaled by
% MRorCT ('mr' or 'ct') Determines whether the template contrast exhibits a CT or MR scan 
%
% Output:
%
% for each scan:
% 
% xCircleFinal Cylinder x-locations after removing the desired cylinders (gnd corrected data)
% yCircleFinal Cylinder y-locations after removing the desired cylinders (gnd corrected data)
% xLocFinal Cylinder x-locations after removing the desired cylinders (cylinder corrected data)
% yLocFinal Cylinder y-locations after removing the desired cylinders (cylinder corrected data)
% xCorrGndFinal Ground-truth x-locations after removing the desired cylinders (gnd corrected data)
% yCorrGndFinal Ground-truth x-locations after removing the desired cylinders (gnd corrected data)
% xGndFinal Ground-truth x-locations after removing the desired cylinders (cylinder corrected data)
% yGndFinal Ground-truth x-locations after removing the desired cylinders (cylinder corrected data)
% radiusFinal The distance from isocenter for all the cylinders after removing cylinders
% totDiffFinal The deviation between the ground truth and cylinder locations after removing cylinders
%
% the combined data:
%
% xLocFinalCombined Cylinder x-locations after removing the desired cylinders (cylinder corrected data)
% yLocFinalCombined Cylinder y-locations after removing the desired cylinders (cylinder corrected data)
% xGndFinalCombined Ground-truth x-locations after removing the desired cylinders (cylinder corrected data)
% yGndFinalCombined Ground-truth x-locations after removing the desired cylinders (cylinder corrected data)
% radiusFinalCombined The distance from isocenter for all the cylinders after removing cylinders
% totDiffFinalCombined The deviation between the ground truth and cylinder locations after removing cylinders
% numRemoved The number of cylinders removed
%
% John Ginn
% Created: 11/16/16
% Modified: 11/17/16

function [xCircleFinal,yCircleFinal,xLocFinal,yLocFinal,...
    xCorrGndFinal,yCorrGndFinal,xGndFinal,yGndFinal,radiusFinal,totDiffFinal,...
    xLocFinalCombined,yLocFinalCombined,xGndFinalCombined,yGndFinalCombined,...
    radiusFinalCombined,totDiffFinalCombined,numRemoved] = ...
    plotOutliers2D(xCircleAll,yCircleAll,finalCircleDataX,finalCircleDataY,deviation2DScans,...
    xCorrGndTruthAll,yCorrGndTruthAll,xGndTruthFinalAll,yGndTruthFinalAll,...
    radiusAllScans,imageDataAll,radiusRemove,rangeDevRemove,sliceShiftExcel,...
    voxelWidthExcel,voxelHeightExcel,radiusSearch1Excel,radiusSearch2Excel,radiusSearch3Excel,scaleFactor,MRorCT)

% initialize the data
xCircleFinal = cell(1,length(xCircleAll));
yCircleFinal = cell(1,length(xCircleAll));
xLocFinal = cell(1,length(xCircleAll));
yLocFinal = cell(1,length(xCircleAll));
xCorrGndFinal = cell(1,length(xCircleAll));
yCorrGndFinal = cell(1,length(xCircleAll));
xGndFinal = cell(1,length(xCircleAll));
yGndFinal= cell(1,length(xCircleAll));
radiusFinal = cell(1,length(xCircleAll));
totDiffFinal = cell(1,length(xCircleAll));

% symbol size
scaleSymbol = ceil(scaleFactor/2);

countFinal = 0;
countHighlight = 0;
% step through all the analyses
for stepScans = 1:length(xCircleAll)
    % reset the data for the scan
    countFinalScan = 0;
    countHighlightScan = 0;
    xLocFinalScanTemp = 0;
    yLocFinalScanTemp = 0;
    xCircleFinalTemp = 0;
    yCircleFinalTemp = 0;
    xCorrGndFinalTemp = 0;
    yCorrGndFinalTemp = 0;
    xGndFinalScanTemp = 0;
    yGndFinalScanTemp = 0;
    radiusFinalScanTemp = 0;
    totDiffFinalScanTemp = 0;
    xDistHighlightScanTemp = 0; 
    yDistHighlightScanTemp = 0;
    xGndHighlightScanTemp = 0;
    yGndHighlightScanTemp = 0;
    radiusHighlightScanTemp = 0;
    totDiffHighlightScanTemp = 0;

    
    % plotting locations for cylinder and ground truth (gnd rot/trans corrected)
    xCylinderPlot = xCircleAll{stepScans};
    yCylinderPlot = yCircleAll{stepScans};
    xGndTruthPlot = xCorrGndTruthAll{stepScans};
    yGndTruthPlot = yCorrGndTruthAll{stepScans};
    % analysis locations for cylinder and ground truth (cylinder rot/trans corrected)
    xCylinder = finalCircleDataX{stepScans};
    yCylinder = finalCircleDataY{stepScans};
    xGndTruth = xGndTruthFinalAll{stepScans};
    yGndTruth = yGndTruthFinalAll{stepScans};
    radius = radiusAllScans{stepScans};
    sliceShift = sliceShiftExcel{stepScans};
    deviation2D = deviation2DScans{stepScans};
    currentImg = imageDataAll{stepScans};
    voxelWidth = voxelWidthExcel{stepScans};
    voxelHeight = voxelHeightExcel{stepScans};
    radiusSearch1 = radiusSearch1Excel{stepScans};
    radiusSearch2 = radiusSearch2Excel{stepScans};
    radiusSearch3 = radiusSearch3Excel{stepScans};
    
    % original information
    origX = 1:length(currentImg(1,:));
    origY = 1:length(currentImg(:,1));
    [origMeshX,origMeshY] = meshgrid(origX,origY);
    % new data
    if (~strcmp(MRorCT,'ct'))
    volUpscaleX = 1/scaleFactor:1/scaleFactor:length(currentImg(1,:));
    volUpscaleY = 1/scaleFactor:1/scaleFactor:length(currentImg(:,1));
    [upMeshX, upMeshY] = meshgrid(volUpscaleX,volUpscaleY);
    % the upscaled data
    SigDataUpscale = interp2(origMeshX,origMeshY,currentImg,...
        upMeshX,upMeshY,'cubic');
    
    upVoxWidth = voxelWidth/scaleFactor;
    upVoxHeight = voxelHeight/scaleFactor;
    else
        SigDataUpscale = currentImg;
        upVoxWidth = voxelWidth;
        upVoxHeight = voxelHeight;
    end
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
    
    % search region boundary plotting
    [circleImg1] = imageCircle(SigDataUpscale,radiusSearch1,upVoxHeight,upVoxWidth,...
        upCenterRow,upCenterCol,sliceShift);
    [circleImg2] = imageCircle(SigDataUpscale,radiusSearch2,upVoxHeight,upVoxWidth,...
        upCenterRow,upCenterCol,sliceShift);
    [circleImg3] = imageCircle(SigDataUpscale,radiusSearch3,upVoxHeight,upVoxWidth,...
        upCenterRow,upCenterCol,sliceShift);
    
    % step through each individual analysis
    for step = 1:length(radius)        
        if (radius(step)>=radiusRemove(1))&&(radius(step)<=radiusRemove(2))...
                &&(deviation2D(step)>=rangeDevRemove(1))&&(deviation2D(step)<=rangeDevRemove(2))
            % reset the image data
            cylinderFitPtsFail = zeros(length(SigDataUpscale(:,1)),length(SigDataUpscale(1,:)));
            cylinderFitPtsPass = zeros(length(SigDataUpscale(:,1)),length(SigDataUpscale(1,:)));
            grndTruthPts = zeros(length(SigDataUpscale(:,1)),length(SigDataUpscale(1,:)));
            % round locations
            if (~strcmp(MRorCT,'ct'))
                data = round([yCylinderPlot(step) xCylinderPlot(step)].*scaleFactor);
                xMarkerLoc = (data(1)-1*scaleSymbol):1:(data(1)+1*scaleSymbol); % make + symbol on image
                yMarkerLoc = (data(2)-1*scaleSymbol):1:(data(2)+1*scaleSymbol);% make + symbol on image
            else
                data = round([yCylinderPlot(step) xCylinderPlot(step)]);
                xMarkerLoc = (data(1)-1*2):1:(data(1)+1*2); % make + symbol on image
                yMarkerLoc = (data(2)-1*2):1:(data(2)+1*2);% make + symbol on image
            end
            % plot the image to determine if you want to remove this datapoint
            % upscale the image for plotting
            %
            % markers based on whether or not in selected range
            cylinderFitPtsFail(xMarkerLoc,yMarkerLoc) = makeShape(length(xMarkerLoc),'x','thick'); % make a x
            
            % store location of ground truth for cylinders
            if (~strcmp(MRorCT,'ct'))
                data = round([yGndTruthPlot(step) xGndTruthPlot(step)].*scaleFactor);
                xMarkerLoc = (data(1)-1*scaleSymbol):1:(data(1)+1*scaleSymbol); % make + symbol on image
                yMarkerLoc = (data(2)-1*scaleSymbol):1:(data(2)+1*scaleSymbol);% make + symbol on image
            else
                data = round([yGndTruthPlot(step) xGndTruthPlot(step)]);
                xMarkerLoc = (data(1)-1):1:(data(1)+1); % make + symbol on image
                yMarkerLoc = (data(2)-1):1:(data(2)+1);% make + symbol on image
            end
            grndTruthPts(xMarkerLoc,data(2)) = 1; % - portion of + marker for ground truth of cylinders
            grndTruthPts(data(1),yMarkerLoc) = 1; % - portion of + marker for ground truth of cylinders
            
            % step through the rest of the data and plot as green squares
            for stepIndImg = 1:length(radius)
                if stepIndImg ~= step
                    if (~strcmp(MRorCT,'ct'))
                        data = round([yCylinderPlot(stepIndImg) xCylinderPlot(stepIndImg)].*scaleFactor);
                        xMarkerLoc = ((data(1)-1*scaleSymbol):1:(data(1)+1*scaleSymbol)); % make + symbol on image
                        yMarkerLoc = ((data(2)-1*scaleSymbol):1:(data(2)+1*scaleSymbol));% make + symbol on image
                    else
                        data = round([yCylinderPlot(stepIndImg) xCylinderPlot(stepIndImg)]);
                        xMarkerLoc = ((data(1)-2):1:(data(1)+2)); % make + symbol on image
                        yMarkerLoc = ((data(2)-2):1:(data(2)+2));% make + symbol on image
                    end
                    % markers based on whether or not in selected range
                    cylinderFitPtsPass(xMarkerLoc,yMarkerLoc) = makeShape(length(xMarkerLoc),'square','thick'); % make a square
                    
                    % store location of ground truth for cylinders
                    if (~strcmp(MRorCT,'ct'))
                        data = round([yGndTruthPlot(stepIndImg) xGndTruthPlot(stepIndImg)].*scaleFactor);
                        xMarkerLoc = (data(1)-1*scaleSymbol):1:(data(1)+1*scaleSymbol); % make + symbol on image
                        yMarkerLoc = (data(2)-1*scaleSymbol):1:(data(2)+1*scaleSymbol);% make + symbol on image
                    else
                        data = round([yGndTruthPlot(stepIndImg) xGndTruthPlot(stepIndImg)]);
                        xMarkerLoc = (data(1)-1):1:(data(1)+1); % make + symbol on image
                        yMarkerLoc = (data(2)-1):1:(data(2)+1);% make + symbol on image
                    end
                    grndTruthPts(xMarkerLoc,data(2)) = 1; % - portion of + marker for ground truth of cylinders
                    grndTruthPts(data(1),yMarkerLoc) = 1; % - portion of + marker for ground truth of cylinders
                end
            end
            % now plot the image showing the location in question
            figure;
            % fit marker plotting
            imshow(SigDataUpscale,[])
            for stepRegions = 1:3
                if stepRegions == 1
                    currentImg = circleImg1;
                elseif stepRegions == 2
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
            set(gca,'XTickMode','manual')
            set(gca,'YTickMode','manual')
            set(gca,'XTick',xAxisLocation)
            set(gca,'YTick',yAxisLocation)
            set(gca,'XTickLabel',xAxisValue)
            set(gca,'YTickLabel',yAxisValue)
            set(gca,'FontSize',16)
            
            if input('Remove this cylinder from the analysis? (y = 1, n = 0)')
                % this location was removed
                % plot as red to show these locations on all the plots
                % store all data in one array
                countHighlight = countHighlight + 1;
                xDistHighlightAll(countHighlight) = xCylinder(step); % index in the image
                yDistHighlightAll(countHighlight) = yCylinder(step);
                xGndHighlightAll(countHighlight) = xGndTruth(step);
                yGndHighlightAll(countHighlight) = yGndTruth(step);
                radiusHighlightAll(countHighlight) = radius(step);
                totDiffHighlightAll(countHighlight) = deviation2D(step);
                % store all data for each scan
                countHighlightScan = countHighlightScan + 1;
                xDistHighlightScanTemp(countHighlightScan) = xCylinder(step); % index in the image
                yDistHighlightScanTemp(countHighlightScan) = yCylinder(step);
                xGndHighlightScanTemp(countHighlightScan) = xGndTruth(step);
                yGndHighlightScanTemp(countHighlightScan) = yGndTruth(step);
                radiusHighlightScanTemp(countHighlightScan) = radius(step);
                totDiffHighlightScanTemp(countHighlightScan) = deviation2D(step);

            else
                % this location was not removed
                % store all data in one array
                countFinal = countFinal + 1;
                xLocFinalCombined(countFinal) = xCylinder(step);
                yLocFinalCombined(countFinal) = yCylinder(step);
                xGndFinalCombined(countFinal) = xGndTruth(step);
                yGndFinalCombined(countFinal) = yGndTruth(step);
                radiusFinalCombined(countFinal) = radius(step);
                totDiffFinalCombined(countFinal) = deviation2D(step);
                % store all data for each scan
                countFinalScan = countFinalScan + 1;
                xCircleFinalTemp(countFinalScan) = xCylinderPlot(step);
                yCircleFinalTemp(countFinalScan) = yCylinderPlot(step);
                xLocFinalScanTemp(countFinalScan) = xCylinder(step);
                yLocFinalScanTemp(countFinalScan) = yCylinder(step);
                xCorrGndFinalTemp(countFinalScan) = xGndTruthPlot(step);
                yCorrGndFinalTemp(countFinalScan) = yGndTruthPlot(step);
                xGndFinalScanTemp(countFinalScan) = xGndTruth(step);
                yGndFinalScanTemp(countFinalScan) = yGndTruth(step);
                radiusFinalScanTemp(countFinalScan) = radius(step);
                totDiffFinalScanTemp(countFinalScan) = deviation2D(step);

            end
            close(gcf) % close the image
        else
            % plot as blue to show these locations on all the plots
            % store all data in one array
            countFinal = countFinal + 1;
            xLocFinalCombined(countFinal) = xCylinder(step);
            yLocFinalCombined(countFinal) = yCylinder(step);
            xGndFinalCombined(countFinal) = xGndTruth(step);
            yGndFinalCombined(countFinal) = yGndTruth(step);
            radiusFinalCombined(countFinal) = radius(step);
            totDiffFinalCombined(countFinal) = deviation2D(step);
            % store all data for each scan
            countFinalScan = countFinalScan + 1;
            xCircleFinalTemp(countFinalScan) = xCylinderPlot(step);
            yCircleFinalTemp(countFinalScan) = yCylinderPlot(step);
            xLocFinalScanTemp(countFinalScan) = xCylinder(step);
            yLocFinalScanTemp(countFinalScan) = yCylinder(step);
            xCorrGndFinalTemp(countFinalScan) = xGndTruthPlot(step);
            yCorrGndFinalTemp(countFinalScan) = yGndTruthPlot(step);
            xGndFinalScanTemp(countFinalScan) = xGndTruth(step);
            yGndFinalScanTemp(countFinalScan) = yGndTruth(step);
            radiusFinalScanTemp(countFinalScan) = radius(step);
            totDiffFinalScanTemp(countFinalScan) = deviation2D(step);
        end
    end
    
    % store the data for each scan
    xCircleFinal{stepScans} = xCircleFinalTemp;
    yCircleFinal{stepScans} = yCircleFinalTemp;
    xLocFinal{stepScans} = xLocFinalScanTemp;
    yLocFinal{stepScans} = yLocFinalScanTemp;
    xCorrGndFinal{stepScans} = xCorrGndFinalTemp;
    yCorrGndFinal{stepScans} = yCorrGndFinalTemp;
    xGndFinal{stepScans} = xGndFinalScanTemp;
    yGndFinal{stepScans} = yGndFinalScanTemp;
    radiusFinal{stepScans} = radiusFinalScanTemp;
    totDiffFinal{stepScans} = totDiffFinalScanTemp;
end






% plot the removed cylinders in red
figure;
subplot(2,2,1)
plot(xLocFinalCombined,totDiffFinalCombined,'.','MarkerSize',8)
if countHighlight > 0
    hold on
    plot(xDistHighlightAll,totDiffHighlightAll,'.r','MarkerSize',8)
end
xlabel('x-image position (index)','FontSize',22)
ylabel('2D Deviation (mm)','FontSize',22)
title('Deviation according to x-position','FontSize',22)
% ylim([0 maxDiff*1.1])
% xlim([0 maxPosX*1.1])
subplot(2,2,2)
plot(yLocFinalCombined,totDiffFinalCombined,'.','MarkerSize',8)
if countHighlight > 0
    hold on
    plot(yDistHighlightAll,totDiffHighlightAll,'.r','MarkerSize',8)
end
xlabel('y-image position (index)','FontSize',22)
ylabel('2D Deviation (mm)','FontSize',22)
title('Deviation according to y-position','FontSize',22)
% ylim([0 maxDiff*1.1])
% xlim([0 maxPosY*1.1])
subplot(2,2,3)
plot(radiusFinalCombined,totDiffFinalCombined,'.','MarkerSize',8)
if countHighlight > 0
    hold on
    plot(radiusHighlightAll,totDiffHighlightAll,'.r','MarkerSize',8)
end
% ylim([0 maxDiff*1.1])
% xlim([0 maxPos*1.1])
xlabel('Distance from Isocenter (index)','FontSize',22)
ylabel('2D Deviation (mm)','FontSize',22)
% title('3D Deviation','FontSize',22)
if exist('totDiffHighlightAll')
    numRemoved = length(totDiffHighlightAll);
else
    numRemoved = 0;
end


end