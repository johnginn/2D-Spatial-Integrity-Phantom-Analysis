% Function to remove the cylinders with the GUI
%
% Input:
% xCircle (index) The found x-location of the cylinders before corrections are applied in each scan (gnd corrected data)
% yCircle (index) The found y-location of the cylinders before corrections are applied in each scan (gnd corrected data)
% xLoc (index) The found x-location of the cylinders after corrections are applied in each scan (cylinder corrected data)
% yLoc (index) The found y-location of the cylinders after corrections are applied in each scan (cylinder corrected data)
% xCorrGnd (index) The ground-truth x-position after corrections for each scan (gnd corrected data)
% yCorrGnd (index) The ground-truth y-position after corrections for each scan (gnd corrected data)
% xGnd (index) The ground-truth x-position before corrections for each scan (cylinder corrected data)
% yGnd (index) The ground-truth y-position before corrections for each scan (cylinder corrected data)
% imageDataAll The image data for all the scans
% sliceShiftExcel The shift of the phantom from isocenter for all scans
% voxelLengthExcel (mm) The voxel length for all scans
% nScans The number of scans
% upscaleFactor The factor to upscale the images by
% MRorCT Whether this is an mr or ct scan 'mr' or 'ct' 
% voxelWidth (mm) The voxel width 
% voxelHeight (mm) The voxel height 
% radiusSearch1 (mm) The radius defining the frist analysis region
% radiusSearch2 (mm) The radius defining the second analysis region
% radiusSearch3 (mm) The radius defining the third analysis region
% deviation2D (mm) Deviation between the ground truth and the found location in 2D
% radius (mm) Distance of the spheres from isocenter
%
% Output:
% numRemoved The number of cylinders removed from the analysis
% xCircleFinal Cylinder x-locations after removing the desired cylinders (gnd corrected data)
% yCircleFinal Cylinder y-locations after removing the desired cylinders (gnd corrected data)
% xLocFinal Cylinder x-locations after removing the desired cylinders (cylinder corrected data)
% yLocFinal Cylinder y-locations after removing the desired cylinders (cylinder corrected data)
% xCorrGndFinal Ground-truth x-locations after removing the desired cylinders (gnd corrected data)
% yCorrGndFinal Ground-truth x-locations after removing the desired cylinders (gnd corrected data)
% xGndFinal Ground-truth x-locations after removing the desired cylinders (cylinder corrected data)
% yGndFinal Ground-truth x-locations after removing the desired cylinders (cylinder corrected data)
% deviation2DFinal (mm) Deviation between the ground truth and the found location in 2D after removing cylinders
% radiusFinal (mm) Distance of the spheres from isocenter after removing cylinders
%
% John Ginn
% Created: 12/8/16
% Modified: 12/13/16


function [numRemoved,xCircleFinal,yCircleFinal,xLocFinal,yLocFinal,...
    xCorrGndFinal,yCorrGndFinal,xGndFinal,yGndFinal,deviation2DFinal,radiusFinal] = ...
    removeCylindersWithGUI_single(xCircle,yCircle,xLoc,yLoc,...
    xCorrGnd,yCorrGnd,xGnd,yGnd,imageDataAll,sliceShiftExcel,voxelLengthExcel,nScans,... 
    upscaleFactor,MRorCT,voxelWidth,voxelHeight,radiusSearch1,radiusSearch2,radiusSearch3,...
    deviation2D,radius)

% cylindersToRemove contains the data in the format [x,y,scan number]
cylindersToRemove = RemoveCylinderGUI_single(xCircle,yCircle,xCorrGnd,yCorrGnd,...
    imageDataAll,sliceShiftExcel,voxelLengthExcel,nScans,... 
    upscaleFactor,MRorCT,voxelWidth,voxelHeight,...
    radiusSearch1,radiusSearch2,radiusSearch3); 

% Extract the analyses by slices
if isempty(cylindersToRemove)
    numRemoved = 0;
else
    numRemoved = length(cylindersToRemove(:,1));
end
xCircleFinal = xCircle;
yCircleFinal = yCircle;
xLocFinal = xLoc;
yLocFinal = yLoc;
xCorrGndFinal = xCorrGnd;
yCorrGndFinal = yCorrGnd;
xGndFinal = xGnd;
yGndFinal = yGnd;
deviation2DFinal = deviation2D;
radiusFinal = radius;
if numRemoved > 0
    for stepCylinders = 1:numRemoved;
        % the current scan number
        thisScan = cylindersToRemove(stepCylinders,3);
        % select the data from each scan
        finalCylinderDataXCorr = xLocFinal{thisScan};
        finalCylinderDataYCorr = yLocFinal{thisScan};
        % calculate the distance to all the cylinders
        distToCylinders = ((finalCylinderDataXCorr.*upscaleFactor - cylindersToRemove(stepCylinders,1)).^2 + ...
            (finalCylinderDataYCorr.*upscaleFactor - cylindersToRemove(stepCylinders,2)).^2);
        % find the closest cylinder
        [minDist, cylinderInd] = min(distToCylinders);
        % extrac out the current Data
        xCircleCurr = xCircleFinal{thisScan};
        yCircleCurr = yCircleFinal{thisScan};
        xLocCurrCurr = xLocFinal{thisScan};
        yLocCurrCurr = yLocFinal{thisScan};
        xCorrGndCurr = xCorrGndFinal{thisScan};
        yCorrGndCurr = yCorrGndFinal{thisScan};
        xGndCurr = xGndFinal{thisScan};
        yGndCurr = yGndFinal{thisScan};
        deviation2DCurr = deviation2DFinal{thisScan};
        radiusCurr = radiusFinal{thisScan};
        % temporary array for truncating the data array
        xCircleTemp = zeros(1,(length(xCircleCurr)-1));
        yCircleTemp = xCircleTemp;
        xLocTemp = xCircleTemp;
        yLocTemp = xCircleTemp;
        xCorrGndTemp = xCircleTemp;
        yCorrGndTemp = xCircleTemp;
        xGndTemp = xCircleTemp;
        yGndTemp = xCircleTemp;
        deviation2DTemp = xCircleTemp;
        radiusTemp = xCircleTemp;
        countTemp = 0;
        for step = 1:length(xCircleCurr)
            if step == cylinderInd
                % skip this cylinder
            else
                countTemp = countTemp + 1;
                xCircleTemp(countTemp) = xCircleCurr(step);
                yCircleTemp(countTemp) = yCircleCurr(step);
                xLocTemp(countTemp) = xLocCurrCurr(step);
                yLocTemp(countTemp) = yLocCurrCurr(step);
                xCorrGndTemp(countTemp) = xCorrGndCurr(step);
                yCorrGndTemp(countTemp) = yCorrGndCurr(step);
                xGndTemp(countTemp) = xGndCurr(step);
                yGndTemp(countTemp) = yGndCurr(step);
                deviation2DTemp(countTemp) = deviation2DCurr(step);
                radiusTemp(countTemp) = radiusCurr(step);
            end
        end
        xCircleFinal{thisScan} = xCircleTemp;
        yCircleFinal{thisScan} = yCircleTemp;
        xLocFinal{thisScan} = xLocTemp;
        yLocFinal{thisScan} = yLocTemp;
        xCorrGndFinal{thisScan} = xCorrGndTemp;
        yCorrGndFinal{thisScan} = yCorrGndTemp;
        xGndFinal{thisScan} = xGndTemp;
        yGndFinal{thisScan} = yGndTemp;
        deviation2DFinal{thisScan} = deviation2DTemp;
        radiusFinal{thisScan} = radiusTemp;
    end
end

% check to see that they are removed
% RemoveCylinderGUI(xCircleFinal,yCircleFinal,xLocFinal,yLocFinal,...
%     imageDataAll,sliceShiftExcel,voxelLengthExcel,nScans,... 
%     upscaleFactor,MRorCT,voxelWidth,voxelHeight,...
%     radiusSearch1,radiusSearch2,radiusSearch3); 

end