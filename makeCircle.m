% This function makes the circle used to step through the images to
% calculate a correlation coefficient
%
% Input:
% voxelHeight (mm/voxel) The height of the pixel
% voxelWidth (mm/voxel) The width of the pixel
% circleDiameter (mm) The diameter of the cylinder
% circleBoundary Whether or not you want the cylinder model to be surrounded
% by a cubic shell of 0's y = 1, n = 0
% MRorCT ('mr' or 'ct') Determines whether the template contrast exhibits a CT or MR scan 
%
% Output:
% circleData (pixels) The image of the circle(1's and 0's) where the cylinder exists
%
% John Ginn
% Created: 11/2/16
% Modified: 11/2/16

function [circleData] = makeCircle(voxelHeight, voxelWidth, circleDiameter,circleBoundary,MRorCT)
circleRadius = circleDiameter/2;

% find necessary range to sample
xRange = round(circleRadius/voxelWidth);
yRange = round(circleRadius/voxelHeight);

% used for storing the cylinder data
countX = 1;
countY = 1;

% temporary cylinder model until it is determine whether or not the cylinder
% needs to be trimmed
circleTemp = zeros(length(-yRange:1:yRange),length(-xRange:1:xRange));
for X = -xRange:1:xRange;
    % reset the count in the matrix
    countY = 1;
    for Y = -yRange:1:yRange;
        % reset the count in the matrix
        % circle equation --> x^2 + y^2 = r^2
        if ((X*voxelWidth)^2 + (Y*voxelHeight)^2 <= circleRadius^2)
            % if inside the cylinder, or on boundary assign a 1
            circleTemp(countY,countX) = 1;
        end
        countY = countY + 1;
    end
    countX = countX + 1;
end

% count the number of ones on the x boundary
countOnesX = 0;
countOnesY = 0;

% dimensions of the circle
yDim = length(circleTemp(:,1,1));
xDim = length(circleTemp(1,:,1));
for X = 1:xDim
    for Y = 1:yDim
        if (X == 1)||(X == xDim)
            if (circleTemp(Y,X) == 1)
                % there is a 1 on one of the x-edges of the volume
                countOnesX = countOnesX + 1;
            end
        end
        if (Y == 1)||(Y == yDim)
            if (circleTemp(Y,X) == 1)
                % there is a 1 on one of the y-edges of the volume
                countOnesY = countOnesY + 1;
            end
        end
    end
end

% determine if there will be a boundary of zeros or not
if circleBoundary == 1 % yes, boundary of zeros
    % need to add zeros to each dimension separately, because the voxel
    % size might not always be symmetric
    if (countOnesX > 0)
        % add a boundary
        xSize = xDim + 2;
        % the location in the new dataset where the old data will be stored
        xArray = (1:xDim) + 1;
    else
        % no boundary needed
        xSize = xDim;
        xArray = (1:xDim);
    end
    if (countOnesY > 0)
        % add a boundary
        ySize = yDim + 2;
        % the location in the new dataset where the old data will be stored
        yArray = (1:yDim) + 1;
    else
        % no boundary needed
        ySize = yDim;
        yArray = (1:yDim);
    end

    % initialize the new data
    circleData = zeros(ySize,xSize);
    circleData(yArray,xArray) = circleTemp;

else % no, no boundary of zeros
    % trim the edges of the data if cylinderRadius/voxelWidth is not an integer
    % value. (avoids using 0's surrounding the cylinder model in the calculation
    % of the correlation coefficient)
    if(countOnesX > 0)
        % don't trim x-dimension
        xData = 1:length(circleTemp(1,:,1));
    else
        % trim the x-dimension
        xData = 2:(length(circleTemp(1,:,1)) - 1);
    end
    if(countOnesY > 0)
        % don't trim y-dimension
        yData = 1:length(circleTemp(:,1,1));
    else
        % trim the y-dimension
        yData = 2:(length(circleTemp(:,1,1)) - 1);
    end
    % store the final data
    xDataFinal = 1:length(xData); % x-dimension may be a different length, start at 1
    yDataFinal = 1:length(yData); % y-dimension may be a different length, start at 1
    circleData(yDataFinal,xDataFinal) = circleTemp(yData,xData);
end

% if you are using a CT scan, invert the contrast
if (strcmp(MRorCT,'ct') == 1)
    for xStep = 1:length(circleData(:,1))
        for yStep = 1:length(circleData(1,:))
            if circleData(xStep,yStep) > 0
                circleData(xStep,yStep) = 0;
            else
                circleData(xStep,yStep) = 1;
            end
        end
    end
end

end