% Script to analyze a single scan of the 2D spatial integrity phantom. This
% script can run the analysis on either a .dat file from the real-time
% imaging protocol or on a dicom image
%
% All the dicom images, or .mat file of the spatial integrity phantom must
% be in a folder without any additional files in that folder
%
% The user must adjust all parameters according to their acquisition parameters
% in order for the analysis code to work appropriately
%
% John Ginn
% Created: 12/1/16
% Modified: 2/7/16

clear all
close all

%% Parameters
% real-time imaging sequence
voxelHeight = 3.5; % (mm) y-dimension (anterior-posterior)
voxelWidth = 3.5; % (mm) x-dimension (medial-lateral)
% (mm) z-dimension (inferior-superior)
voxelLength = 5; 
sliceShift = 0;

%% Data information
realTime = 0; % y = 1, n = 0 Whether or not you are using images from the "real time" pulse sequence
% realTime requires the .dat file from the scanner where as the 
obtainFilenames = 1; % y = 1, n = 0 Whether or not you need to obtain the filenames for the real time imaging

% fileName = filenameSorted{analysisNumber}; % data file
numOfFrames = 40; % number of frames acquired in the .dat file if the real-time imaging sequence was used

%% Options
plotData = 1; % y = 1, n = 0 Whether or not to plot the data
plotOutsideRegion = 0; % y = 1, n = 0 Whether or not to plot the cylinders passing or failing radiusThershol3 outside radius 3
corrThreshold = 0.65; % the minimum correlation threshold
radiusSearch1 = 50;  % (mm) distance for searching for cylinders within 100 mm of isocenter
radiusThreshold1 = 1; % (mm) the passing criteria for cylinders within 100 mm of isocenter
radiusSearch2 = 100;  % (mm) distance for searching for cylinders within 175 mm of isocenter
radiusThreshold2 = 1; % (mm) the passing criteria for cylinders within 175 mm of isocenter
radiusSearch3 = 175;  % (mm) distance for searching for cylinders within 100 mm of isocenter
radiusThreshold3 = 2; % (mm) the passing criteria for cylinders within 100 mm of isocenter
upscaleFactor = 7;
MRorCT = 'mr';

% do you need to obtain the center slice, center column and center row location?
obtainLocations = 1; % y = 1,n = 0
% WARNING: this will overwrite images saved in this direcory with the same name
saveExcelFiles = 1; % y = 1, n = 0 save the data to an excel file
saveImages = 0; % y = 1, n = 0 save the analysis images at a high resolution 
% saveImageDirectory = strcat([pwd,'\AnalysisImages']);
saveImageDirectory = pwd;
analysisType = 'upscale';
% 'upscale' if you want to upscale the data for the analysis
% use the 'weighted' method that calculates
% the location of the cylinders based on the weighted sum of the signal
% Or 'non-upscale' locations based off the non-upscaled correlation coefficient
% Note: upscaling can greatly increase the computational time
saveDicom = 0; % whether or not to create and save the dicom images (raw data images, not analysis)

if saveImages == 1
    plotData = 1;
end
%% Other Options and Parameters Rarely Changed 
% parameters
nScanDirectories = 1; % the number of different scan directories containing data
searchDist = 4; % (mm) distance in each direction to search for the cylinder (orig corr coeff)
searchDistWeight = 4; % (mm) distance in each direction to use to calulate weighted sum location
binSize = 20; % (mm) The size of the bins for plotting the deviation
circleDiameter = 13;

% options
circleBoundary = 1; % y = 1,n = 0 make cylinder model have a boundary of zeros around it

if strcmp(analysisType,'upscale'); 
    % y = 1, n = 0 construct an upscaled dataset and compute the correlation coefficent to find the location of the cylinder
    upscale = 1;
else
    upscale = 0;
end
%
% change if you already know what these values are,
% otherwise they will be overwritten when you run the code
if obtainLocations == 0
    tenthCol = 12;
    tenthRow = 11;
end

% phantom information
cylinderDiameter = 8; % (mm)
distBtwnCylinders = 14.57; % (mm) The distance between the cylinders in each direction
xIndPerSph = distBtwnCylinders/voxelWidth; % number of pixels in x-direction between the cylinders
yIndPerSph = distBtwnCylinders/voxelHeight; % number of pixels in y-direction between the cylinders
nRows = 20; % there are 20 rows of cylinders
nTotData = 397; % total number of cylinders

%% Extract the data
nScans = 0;
nScansPerFolder = zeros(1,nScanDirectories);
filenameSortedFolder = cell(1,nScanDirectories);
fileDirectoryFolder = cell(1,nScanDirectories);
for cycleDatasets = 1:nScanDirectories
    % "real time" data information
    if realTime == 1
        [filenameSorted,sortedTime,fileDirectory] = readFilenames();
        nScans = nScans + length(sortedTime);
        % number of scans in each folder
        nScansPerFolder(cycleDatasets) = length(sortedTime);
        filenameSortedFolder{cycleDatasets} = filenameSorted;
        fileDirectoryFolder{cycleDatasets} = fileDirectory;
    else
        currentDir = pwd;
        disp('Please select folder containing the data files')
        fileDirectory = uigetdir;
        cd(fileDirectory)
        filesStructure = dir;
        cd(currentDir);
        fileDirectoryFolder{cycleDatasets} = fileDirectory;
        filenameSortedFolder{cycleDatasets} = fileDirectory;
        count = 0; % then number of data files
        % sort through files and use only .dat files
        for step = 1:(length(filesStructure) - 2)
            % +2 to ignore the first two lines '.' and '..' respectively
            currentFile = filesStructure(step + 2).name;
            if strcmp('.dat',currentFile((length(currentFile)-3):length(currentFile)))||...
                    strcmp('.IMA',currentFile((length(currentFile)-3):length(currentFile)))
                count = count + 1;
                fileNames{count} = currentFile;
            end
        end
        nScans = 1;
    end
end
missingLocationsFolder = cell(1,nScanDirectories);
radiusAllScans = cell(1,nScans);
deviation2DAllScans = cell(1,nScans);
xDeviationAllScans =  cell(1,nScans);
yDeviationAllScans =  cell(1,nScans);
optCorrAllScans =  cell(1,nScans);

corrThresholdExcel = cell(1,nScans);
radiusThreshold1Excel = cell(1,nScans);
radiusThreshold2Excel = cell(1,nScans);
radiusThreshold3Excel = cell(1,nScans);
percentPassAllExcel = cell(1,nScans);
radiusSearch1Excel = cell(1,nScans);
radiusSearch2Excel = cell(1,nScans);
radiusSearch3Excel = cell(1,nScans);
percentPassXY1Excel = cell(1,nScans);
percentPassXY2Excel = cell(1,nScans);
percentPassXY3Excel = cell(1,nScans);
avgDevRadius1XYExcel = cell(1,nScans);
avgDevRadius2XYExcel = cell(1,nScans);
avgDevRadius3XYExcel = cell(1,nScans);
avgDevAllXYExcel = cell(1,nScans);
maxDevRadius1XYExcel = cell(1,nScans);
maxDevRadius2XYExcel = cell(1,nScans);
maxDevRadius3XYExcel = cell(1,nScans);
maxDevAllExcel = cell(1,nScans);
sliceShiftExcel = cell(1,nScans);
voxelWidthExcel = cell(1,nScans);
voxelHeightExcel = cell(1,nScans);
voxelLengthExcel = cell(1,nScans);
fileNameExcel = cell(1,nScans);
finalCircleDataXAll = cell(1,nScans);
finalCircleDataYAll = cell(1,nScans);
xGndTruthFinalAll = cell(1,nScans);
yGndTruthFinalAll = cell(1,nScans);
imageDataAll = cell(1,nScans);
xCenterAll = cell(1,nScans);
yCenterAll = cell(1,nScans);
xCircleAll = cell(1,nScans);
yCircleAll = cell(1,nScans);
xCorrGndTruthAll = cell(1,nScans);
yCorrGndTruthAll = cell(1,nScans);

countScanDirectories = 1; % for stepping through different folders of data
stepFiles = 1; % count for stepping through the files in each folder

for cycleScans = 1:nScans    
    if nScans > 1;
       disp(' ')
       disp(strcat(['Analysis Count: ' ,num2str(cycleScans)]))
       disp(' ')
    end
    % determine direction of shift for saving plots
    if sign(sliceShift) == 1
        posOrNeg = 'Pos';
    elseif sign(sliceShift) == 0
        posOrNeg = 'Iso';
    else
        posOrNeg = 'Neg';
    end
    % image acquisition using the real time images
    currentDir = pwd;
    if realTime == 1
        fileName = filenameSorted{stepFiles}; % data file
        imageData = extractVolume(fileName,fileDirectory,numOfFrames,saveDicom);
    else
        linesToSkip = 3;
        [fileNames, fileData, fileMap, fileInfo] = loadImages(fileDirectory,currentDir,linesToSkip);
        % Construct the volumetric signal
        volData = zeros(size(fileData{1},1),size(fileData{1},2),length(fileData));
        for step = 1:length(fileData);
            volData(:,:,step) = fileData{step};
        end
        %
        % Plot the volume extracted from images or simulation and obtain locations
        % of the center slice, row and column if necessary
        % store variables for GUI
        guiData{1} = volData;
        guiData{2} = fileMap;
        % display the GUI
        UserInputGUI2(guiData);
        analysisSlice = input('Enter image number for analysis (if single image in folder enter 1): ');
        imageData = volData(:,:,analysisSlice);
        fileName = 'not real-time imaging';
    end
    % If on the first scan, find the missing spheres
    if cycleScans == 1;
        for cycleDatasets = 1:nScanDirectories
        % obtain the appropriate filenames and folder directory
        filenameSortedTemp = filenameSortedFolder{cycleDatasets};
        fileDirectoryTemp = fileDirectoryFolder{cycleDatasets};
        if sum(class(filenameSortedTemp)=='cell') == 4
            fileNameTemp = filenameSortedTemp{cycleScans}; % data file
        else
            fileNameTemp = filenameSortedTemp;
        end
        if realTime == 1
            imageDataTemp = extractVolume(fileNameTemp,fileDirectoryTemp,numOfFrames,saveDicom);
        else
            imageDataTemp = imageData;
        end
            if obtainLocations == 1
                figure;
                imshow(imageDataTemp,[])
                disp(' ')
                disp(' ')
                disp(' ')
                disp('The first cylinder is defined in the upper left-hand corner of the phantom')
                tenthCol  = input('Enter 10th column index, counting left to right (center of the cylinder) [x]: '); % location of center column of cylinders
                tenthRow = input('Enter 10th row index, counting top to bottom (center of the cylinder) [y]: '); % index of the center row
            end
            xDim = length(imageDataTemp(1,:));
            yDim = length(imageDataTemp(:,1));
            [X_ns,Y_ns]=meshgrid((1:xDim)*voxelWidth,(1:yDim)*voxelHeight);
            [Xq,Yq]=meshgrid((1:xDim*upscaleFactor)*voxelWidth/upscaleFactor,(1:yDim*upscaleFactor)*voxelHeight/upscaleFactor);
            upscaleImg = interp2(X_ns,Y_ns,imageDataTemp,Xq,Yq,'cubic');
            % determine the cylinders to skip depending on the orientation of the phantom
            % use the upscaled image to improve accuracy
            % start at zero so first location is first column and row
            firstCol = tenthCol - 9*xIndPerSph;
            firstRow = tenthRow - 9*yIndPerSph;
            possibleXPositionsUp = ((0:19).*xIndPerSph + firstCol).*upscaleFactor;
            possibleYPositionsUp = ((0:19).*yIndPerSph + firstRow).*upscaleFactor;
            figure;
            imshow(upscaleImg,[])
            title('Select Missing Cylinders','FontSize',20)
            disp('Select 1st missing cylinder location')
            [xPos1,yPos1] = ginput(1);
            % find the row and coumn number
            [xPos1min, missingCol1] = min(abs(possibleXPositionsUp - xPos1));
            [yPos1min, missingRow1] = min(abs(possibleYPositionsUp - yPos1));
            disp(strcat(['1st missing cylinder row and column number [row, column]: [',...
                num2str(missingRow1),',',num2str(missingCol1),']']))
            disp('Select 2nd missing cylinder location')
            [xPos2,yPos2] = ginput(1);
            % find the row and coumn number
            [xPos2min, missingCol2] = min(abs(possibleXPositionsUp - xPos2));
            [yPos2min, missingRow2] = min(abs(possibleYPositionsUp - yPos2));
            disp(strcat(['2nd missing cylinder row and column number [row, column]: [',...
                num2str(missingRow2),',',num2str(missingCol2),']']))
            disp('Select 3rd missing cylinder location')
            [xPos3,yPos3] = ginput(1);
            % find the row and coumn number
            [xPos3min, missingCol3] = min(abs(possibleXPositionsUp - xPos3));
            [yPos3min, missingRow3] = min(abs(possibleYPositionsUp - yPos3));
            disp(strcat(['3rd missing cylinder row and column number [row, column]: [',...
                num2str(missingRow3),',',num2str(missingCol3),']']))
            close(gcf)
            missingLocations = [missingRow1,missingCol1;missingRow2,missingCol2;missingRow3,missingCol3];
            missingLocationsFolder{cycleDatasets} = missingLocations;
        end
    end
    % Create Circle template and find cylinders
    [circleData] = makeCircle(voxelHeight, voxelWidth, circleDiameter,circleBoundary,MRorCT);
%     upscale the image
    xDim = length(imageData(1,:));
    yDim = length(imageData(:,1));
    [X_ns,Y_ns]=meshgrid((1:xDim)*voxelWidth,(1:yDim)*voxelHeight);
    [Xq,Yq]=meshgrid((1:xDim*upscaleFactor)*voxelWidth/upscaleFactor,(1:yDim*upscaleFactor)*voxelHeight/upscaleFactor);
    upscaleImg = interp2(X_ns,Y_ns,imageData,Xq,Yq,'cubic');
    
    % update the missing locations depending on the day the  phantom was
    % scanned
    [coordcylinderData, coordcylinderPlot, optCorrelation, cylinderGroundTruth,...
        weightCoord, weightCoordPlot, coordcylinderDataUp,coordcylinderPlotUp, optCorrelationUp,...
        imgSearchRegion,imgAllGndTruth,imgWeightSum] = PhantomScan2D(imageData,upscaleImg, circleData, voxelHeight,voxelWidth,...
        tenthCol,tenthRow,nRows,searchDist,searchDistWeight,circleBoundary,MRorCT,upscaleFactor,missingLocations,distBtwnCylinders);
        
    % Analyze the data
    
    % sort the data into individual slices (required to fit a plane to the data
    % to obtain the rotation of the phantom setup)
    if strcmp(analysisType,'upscale') == 1
        % upscale analysis has occurred
        corrData = optCorrelationUp;
    elseif strcmp(analysisType,'non-upscale') == 1
        % non-upscaled correlation coefficient results
        corrData = optCorrelation;
    else
        % the weighted sum method
        corrData = optCorrelation;
    end
    
    % not meeting correlation threshold
    % initialize arrays
    cylImg = zeros(length(imageData(:,1)),length(imageData(1,:)));
    cylinderFitPts = cylImg; % store just the location of the cylinders
    grndTruthPts = cylImg; % the location of the "ground truth" points
    cylImgX = length(circleData(1,:)); % x-dimension of the cylinder
    cylImgY = length(circleData(:,1)); % y-dimension of the cylinder
    % count the data that actually passes the correlation coeff. threshold
    countPass = 0;
    % count data in each slice that passes
    countSlicePass = 0;
    stepCount = 1; % total step counter
    countFail = 0; % data that fail
    clear xCircle yCircle xGndTruthFinal yGndTruthFinal xCylinderFail yCylinderFail...
    xGndTruthFail yGndTruthFail
    % step through the slices
    for stepCircle = 1:nTotData
        if corrData{stepCount} < corrThreshold
            cylinderLocData = coordcylinderPlot{stepCount};
            GroundTruthFinal = cylinderGroundTruth{stepCount};
            % if correlation is less than the minimum correlation threshold,
            % do not plot the cylinder
            countFail = countFail + 1;
            xCylinderFail(countFail) = cylinderLocData(2);
            yCylinderFail(countFail) = cylinderLocData(1);
            xGndTruthFail(countFail) = GroundTruthFinal(2);
            yGndTruthFail(countFail) = GroundTruthFinal(1);
        else
            % count the data that actually passes the correlation coeff.
            % threshold
            countPass = countPass + 1;
            % correlation threshold met
            % store data for plotting
            if mod(cylImgX,2) == 0; % x-dimension of cylinder is even
                cylRangeX = cylImgX/2;
            else % x-dimension of cylinder is odd
                cylRangeX = (cylImgX-1)/2;
            end
            if mod(cylImgY,2) == 0; % y-dimension of cylinder is even
                cylRangeY = cylImgY/2;
            else % y-dimension of cylinder is odd
                cylRangeY = (cylImgY-1)/2;
            end
            data = coordcylinderPlot{stepCount};
            xLoc = (data(2)-cylRangeX):1:(data(2)+cylRangeX);
            yLoc = (data(1)-cylRangeY):1:(data(1)+cylRangeY);
            xMarkerLoc = (data(1)-1):1:(data(1)+1); % make + symbol on image
            yMarkerLoc = (data(2)-1):1:(data(2)+1);% make + symbol on image
            cylImg(yLoc,xLoc) = circleData; % store the cylinder volume
            cylinderFitPts(xMarkerLoc,data(2)) = 1; % - portion of + marker for fit of cylinders
            cylinderFitPts(data(1),yMarkerLoc) = 1; % | portion of + marker for fit of cylinders
            % store final data without cylinders that did not meet correlation
            % requirement
            % data for cylinder agreement calculation
            if stepCount == 1
                optCorrelationFinal = 0; % clear the optimal correlation array
            end
            if strcmp(analysisType,'upscale') == 1 % upscale analysis has occurred
                cylinderLocData = coordcylinderDataUp{stepCount};
                optCorrelationFinal(countPass) = optCorrelationUp{stepCount};
            elseif strcmp(analysisType,'non-upscale') == 1
                % non-upscaled correlation coefficient results
                cylinderLocData = coordcylinderPlot{stepCount};
                optCorrelationFinal(countPass) = corrData{stepCount};
            else
                % the weighted sum results (uses correlation coefficient
                % from non-upscaled data to remove cylinders with a correlation
                % coefficient less than the specified threshold)
                cylinderLocData = weightCoordPlot{stepCount};
                optCorrelationFinal(countPass) = corrData{stepCount};
            end
            xCircle(countPass) = cylinderLocData(2); % cylinder location x-component
            yCircle(countPass) = cylinderLocData(1); % cylinder location y-component
            % store location of ground truth for cylinders
            data = cylinderGroundTruth{stepCount};
            xMarkerLoc = round(data(1)-1):1:(data(1)+1); % make + symbol on image
            yMarkerLoc = round(data(2)-1):1:(data(2)+1);% make + symbol on image
            grndTruthPts(xMarkerLoc,round(data(2))) = 1; % - portion of + marker for ground truth of cylinders
            grndTruthPts(round(data(1)),yMarkerLoc) = 1; % - portion of + marker for ground truth of cylinders
            % ground truth points that did not meet correlation requirement
            GroundTruthFinal = cylinderGroundTruth{stepCount};
            xGndTruthFinal(countPass) = GroundTruthFinal(2);
            yGndTruthFinal(countPass) = GroundTruthFinal(1);
        end
        stepCount = stepCount + 1; % add one step
    end
    percentPassAll = countPass/nTotData*100;
    
    % Procrustes Method
    % apply the series of rotation and translation corrections
    plotPlaneFitCard = 0;
    procylinderBefore = [xCircle',yCircle'];
    proGround = [xGndTruthFinal',yGndTruthFinal'];
    [d, procylinderAfter,transform] = procrustes(proGround,procylinderBefore,'scaling',false);
    proX = procylinderAfter(:,1);
    proY = procylinderAfter(:,2);
    finalCircleDataX = proX';
    finalCircleDataY = proY';
    
    
    
    % Plotting
    % calculate the difference between the found location and the ground truth
    % locations
    nDataFinal = length(finalCircleDataX);
    deviation2D = zeros(1,nDataFinal);
    xDeviation = deviation2D;
    yDeviation = deviation2D;
    % calculate distance from center of the image
    radius = deviation2D;
    % check if dimensions are even or odd. If odd, just divide by 2 (center
    % position is an integer value), otherwise add 0.5
    if mod(length(imageData(1,:)),2) ~= 0
        xCenter = length(imageData(1,:))/2;
    else
        xCenter = length(imageData(1,:))/2 + 0.5;
    end
    if mod(length(imageData(1,:)),2) ~= 0
        yCenter = length(imageData(:,1))/2;
    else
        yCenter = length(imageData(:,1))/2 + 0.5;
    end
    for step = 1:length(xCircle)
        radius(step) = sqrt( (voxelWidth*(xCenter -  finalCircleDataX(step)))^2 + ...
            (voxelHeight*(yCenter -  finalCircleDataY(step)))^2 + sliceShift.^2);
        deviation2D(step) = sqrt( (voxelWidth*(finalCircleDataX(step) -  xGndTruthFinal(step)))^2 + ...
            (voxelHeight*(finalCircleDataY(step) -  yGndTruthFinal(step)))^2);
        xDeviation(step) =  sqrt( (voxelWidth*(finalCircleDataX(step) -  xGndTruthFinal(step)))^2);
        yDeviation(step) =  sqrt( (voxelHeight*(finalCircleDataY(step) -  yGndTruthFinal(step)))^2);
    end
    
    
    % Shift the ground truth data to center at center of cylinders
    % extract out the rotation and translation information
    proT = transform.T; % rotation matrix from the transform
    proC = transform.c;
    
    
    % apply the reverse of procrustes --> Y = (Z - C)T^-1
    groundBeforeCorr = [xGndTruthFinal',yGndTruthFinal'];
    proCtemp = proC(1,:);
    xProC = proCtemp(1,1).*ones(length(xGndTruthFinal),1);
    yProC = proCtemp(1,2).*ones(length(xGndTruthFinal),1);
    proCsmall = [xProC,yProC];
    gndCorrected = (groundBeforeCorr - proCsmall)/(proT); % this does the inverse

    avgDiffGroundCorr = 0;
    avgDiffCylCorr = 0;
    for step = 1:length(gndCorrected(:,1))
        xCorrGndTruth(step) = gndCorrected(step,1);
        yCorrGndTruth(step) = gndCorrected(step,2);
        % after correcting the ground truth data instead
        avgDiffGroundCorr = avgDiffGroundCorr + sqrt((xCorrGndTruth(step)-xCircle(step)).^2 + ...
            (yCorrGndTruth(step)-yCircle(step)).^2);
        avgDiffCylCorr = avgDiffCylCorr + sqrt((finalCircleDataX(step)-xGndTruthFinal(step)).^2 + ...
            (finalCircleDataY(step)-yGndTruthFinal(step)).^2);
    end
    
    % plot passing and failing locations on top of the image
    [compareDeviationXY,avgDevAllXY,percentPassXY1,percentPassXY2,percentPassXY3,...
        avgDevRadius1XY,avgDevRadius2XY,avgDevRadius3XY,...
        maxDevRadius1XY,maxDevRadius2XY,maxDevRadius3XY]  =...
        compViewRay2D(finalCircleDataX,finalCircleDataY,xGndTruthFinal,yGndTruthFinal,...
        voxelHeight,voxelWidth,imageData,radiusThreshold1,radiusThreshold2,radiusThreshold3,...
        xCenter,yCenter,radiusSearch1,radiusSearch2,radiusSearch3,upscaleFactor,sliceShift,...
        xCircle,yCircle,xCorrGndTruth,yCorrGndTruth,plotData,plotOutsideRegion,MRorCT);
    maxDevAll = max(deviation2D);
    if saveImages == 1;
       cd(saveImageDirectory)
       fig = gcf;
       fig.PaperPositionMode = 'auto';
       plotTitle = strcat(['compareViewRay',posOrNeg,num2str(abs(sliceShift)),...
           'SliceThick',num2str(round(voxelLength))]);
       print(plotTitle,'-dtiff','-r300')
       cd(currentDir)
       close(fig)
    end

    
    % plot deviation with respect to distance from center of the image
    if plotData == 1
        figure;
        plot(radius,deviation2D,'.')
        xlabel('Distance from Isosenter (mm)','FontSize',16)
        ylabel('2D Deviation (mm)','FontSize',16)
        title('2D Deviation from Gound Truth','FontSize',16)
    end
    if saveImages == 1;
        cd(saveImageDirectory)
        fig = gcf;
        fig.PaperPositionMode = 'auto';
        plotTitle = strcat(['IndDevPlotShift',posOrNeg,num2str(abs(sliceShift)),...
            'SliceThick',num2str(round(voxelLength))]);
        print(plotTitle,'-dtiff','-r300')
        cd(currentDir)
        close(fig)
    end
    
    % plot deviation in bins
    if plotData == 1
        plotDeviation(radius,deviation2D,binSize)
    end
    if saveImages == 1;
       cd(saveImageDirectory)
       fig = gcf;
       fig.PaperPositionMode = 'auto';
       plotTitle = strcat(['DevPlotShift',posOrNeg,num2str(abs(sliceShift)),...
           'SliceThick',num2str(round(voxelLength))]);
       print(plotTitle,'-dtiff','-r300')
       cd(currentDir)
       close(fig)
    end
    if nScans == 1;
        disp(' ')
        disp(' ')
        disp(' ')
        disp(strcat(['Percent of cylinders passing correlation threshold of ',num2str(corrThreshold),...
            ': ',num2str(percentPassAll)]))
        disp(strcat(['Cylinders w/in ',num2str(radiusSearch1),' (mm) of image center,'...
            ' less than ',num2str(radiusThreshold1),' (mm) deviation: ',num2str(percentPassXY1),'%']))
        disp(strcat(['Cylinders w/in ',num2str(radiusSearch2),' (mm) of image center,'...
            ' less than ',num2str(radiusThreshold2),' (mm) deviation: ',num2str(percentPassXY2),'%']))
        disp(strcat(['Cylinders w/in ',num2str(radiusSearch3),' (mm) of image center,'...
            ' less than ',num2str(radiusThreshold3),' (mm) deviation: ',num2str(percentPassXY3),'%']))
        disp(strcat(['Cylinders w/in ',num2str(radiusSearch1),' (mm) of image center,'...
            ' avg 2D deviation: ',num2str(avgDevRadius1XY),' (mm)']))
        disp(strcat(['Cylinders w/in ',num2str(radiusSearch2),' (mm) of image center,'...
            ' avg 2D deviation: ',num2str(avgDevRadius2XY),' (mm)']))
        disp(strcat(['Cylinders w/in ',num2str(radiusSearch3),' (mm) of image center,'...
            ' avg 2D deviation: ',num2str(avgDevRadius3XY),' (mm)']))
        disp(strcat(['Avg 2D deviation: ',num2str(avgDevAllXY),' (mm) for all cylinders']))
        disp(strcat(['Max 2D deviation: ',num2str(maxDevRadius1XY),' (mm) for cylinders w/in '...
            ,num2str(radiusSearch1),' (mm) of image center,']))
        disp(strcat(['Max 2D deviation: ',num2str(maxDevRadius2XY),' (mm) for cylinders w/in '...
            ,num2str(radiusSearch2),' (mm) of image center,']))
        disp(strcat(['Max 2D deviation: ',num2str(maxDevRadius3XY),' (mm) for cylinders w/in '...
            ,num2str(radiusSearch3),' (mm) of image center,']))
        disp(strcat(['Max 2D deviation: ',num2str(maxDevAll),' (mm) for all cylinders']))
    end
    % store the different radius locations and deviations
    radiusAllScans{cycleScans} = radius;
    deviation2DAllScans{cycleScans} = deviation2D;
    xDeviationAllScans{cycleScans} =  xDeviation;
    yDeviationAllScans{cycleScans} =  yDeviation;
    optCorrAllScans{cycleScans} = optCorrelationFinal;
    finalCircleDataXAll{cycleScans} = finalCircleDataX;
    finalCircleDataYAll{cycleScans} = finalCircleDataY;
    xGndTruthFinalAll{cycleScans} = xGndTruthFinal;
    yGndTruthFinalAll{cycleScans} = yGndTruthFinal;
    imageDataAll{cycleScans} = imageData;
    xCenterAll{cycleScans} = xCenter;
    yCenterAll{cycleScans} = yCenter;
    xCircleAll{cycleScans} = xCircle;
    yCircleAll{cycleScans} = yCircle;
    xCorrGndTruthAll{cycleScans} = xCorrGndTruth;
    yCorrGndTruthAll{cycleScans} = yCorrGndTruth;

    % store the data for saving the excel files
    corrThresholdExcel{cycleScans} = corrThreshold;
    radiusThreshold1Excel{cycleScans} = radiusThreshold1;
    radiusThreshold2Excel{cycleScans} = radiusThreshold2;
    radiusThreshold3Excel{cycleScans} = radiusThreshold3;
    percentPassAllExcel{cycleScans} = percentPassAll;
    radiusSearch1Excel{cycleScans} = radiusSearch1;
    radiusSearch2Excel{cycleScans} = radiusSearch2;
    radiusSearch3Excel{cycleScans} = radiusSearch3;
    percentPassXY1Excel{cycleScans} = percentPassXY1;
    percentPassXY2Excel{cycleScans} = percentPassXY2;
    percentPassXY3Excel{cycleScans} = percentPassXY3;
    avgDevRadius1XYExcel{cycleScans} = avgDevRadius1XY;
    avgDevRadius2XYExcel{cycleScans} = avgDevRadius2XY;
    avgDevRadius3XYExcel{cycleScans} = avgDevRadius3XY;
    avgDevAllXYExcel{cycleScans} = avgDevAllXY;
    maxDevRadius1XYExcel{cycleScans} = maxDevRadius1XY;
    maxDevRadius2XYExcel{cycleScans} = maxDevRadius2XY;
    maxDevRadius3XYExcel{cycleScans} = maxDevRadius3XY;
    maxDevAllExcel{cycleScans} = maxDevAll;
    sliceShiftExcel{cycleScans} = sliceShift;
    voxelWidthExcel{cycleScans} = voxelWidth;
    voxelHeightExcel{cycleScans} = voxelHeight;
    voxelLengthExcel{cycleScans} = voxelLength;
    fileNameExcel{cycleScans} = fileName;
    stepFiles = stepFiles + 1;
end


%% Check the data for outliers and remove if necessary

xCircleFinal = xCircleAll;
yCircleFinal = yCircleAll;
xLocFinal = finalCircleDataXAll;
yLocFinal = finalCircleDataYAll;
xCorrGndFinal = xCorrGndTruthAll;
yCorrGndFinal = yCorrGndTruthAll;
xGndFinal = xGndTruthFinalAll;
yGndFinal = yGndTruthFinalAll;

[numRemoved,xCircleFinal,yCircleFinal,xLocFinal,yLocFinal,...
    xCorrGndFinal,yCorrGndFinal,xGndFinal,yGndFinal,deviation2DAllScans,radiusAllScans] = ...
    removeCylindersWithGUI_single(xCircleFinal,yCircleFinal,xLocFinal,yLocFinal,...
    xCorrGndFinal,yCorrGndFinal,xGndFinal,yGndFinal,imageDataAll,sliceShiftExcel,voxelLengthExcel,nScans,... 
    upscaleFactor,MRorCT,voxelWidth,voxelHeight,radiusSearch1,radiusSearch2,radiusSearch3,deviation2DAllScans,radiusAllScans);



% combine the data
radiusAll = 0;
deviation2DAll = 0;
for step = 1:nScans
    prevPos = length(radiusAll);
    currData = radiusAllScans{step};
    currentPos = prevPos + length(currData);
    if step == 1
        % at first iteration, length of data is the length of the current
        % data. prevPos adds 1 so 1 must be removed at the end
        dataLoc = prevPos:(currentPos - 1);
    else
        % add one at the beginning to avoid overwriting the data
        dataLoc = (prevPos + 1):currentPos;
    end
    radiusAll(dataLoc) = radiusAllScans{step};
    deviation2DAll(dataLoc) = deviation2DAllScans{step};
end

figure;
plot(radiusAll,deviation2DAll,'.')
xlabel('Distance from Isocenter (mm)','FontSize',16)
ylabel('2D Deviation (mm)','FontSize',16)
title('2D Deviation from Ground Truth, All Analyses','FontSize',22)

countRemoveIteration = 0;
totRemoved = numRemoved;
if input('Do you need to check for cylinders in artifact regions? (y = 1, n = 0):')
    outliersRemoved = 0;
    while outliersRemoved == 0
        radiusRemove(1) = input('Input min radius from isocenter to be checked:');
        radiusRemove(2) = input('Input max radius from isocenter to be checked:');
        rangeDevRemove(1) = input('Input min deviation to be checked:');
        rangeDevRemove(2) = input('Input max deviation to be checked:');
        countRemoveIteration = countRemoveIteration + 1;
        % close the deviation plot
        close(gcf)
        % don't go through all the cylinders again if not on the first
        % round of removing the cylinders
        if countRemoveIteration == 1;
            [xCircleFinal,yCircleFinal,xLocFinal,yLocFinal,...
                xCorrGndFinal,yCorrGndFinal,xGndFinal,yGndFinal,radiusFinal,totDiffFinal,...
                xLocFinalCombined,yLocFinalCombined,xGndFinalCombined,yGndFinalCombined,...
                radiusFinalCombined,totDiffFinalCombined,numRemoved] = ...
                plotOutliers2D(xCircleAll,yCircleAll,finalCircleDataXAll,finalCircleDataYAll,deviation2DAllScans,...
                xCorrGndTruthAll,yCorrGndTruthAll,xGndTruthFinalAll,yGndTruthFinalAll,...
                radiusAllScans,imageDataAll,radiusRemove,rangeDevRemove,sliceShiftExcel,...
                voxelWidthExcel,voxelHeightExcel,radiusSearch1Excel,radiusSearch2Excel,radiusSearch3Excel,upscaleFactor,MRorCT);
            totRemoved = numRemoved + totRemoved;
        else
            [xCircleFinal,yCircleFinal,xLocFinal,yLocFinal,...
                xCorrGndFinal,yCorrGndFinal,xGndFinal,yGndFinal,radiusFinal,totDiffFinal,...
                xLocFinalCombined,yLocFinalCombined,xGndFinalCombined,yGndFinalCombined,...
                radiusFinalCombined,totDiffFinalCombined,numRemoved] = ...
                plotOutliers2D(xCircleFinal,yCircleFinal,xLocFinal,yLocFinal,totDiffFinal,...
                xCorrGndFinal,yCorrGndFinal,xGndFinal,yGndFinal,...
                radiusFinal,imageDataAll,radiusRemove,rangeDevRemove,sliceShiftExcel,...
                voxelWidthExcel,voxelHeightExcel,radiusSearch1Excel,radiusSearch2Excel,radiusSearch3Excel,upscaleFactor,MRorCT);
            totRemoved = numRemoved + totRemoved;
        end
        %  plot the deviation plot again
        figure;
        plot(radiusFinalCombined,totDiffFinalCombined,'.')
        xlabel('Distance from Isocenter (mm)','FontSize',16)
        ylabel('2D Deviation (mm)','FontSize',16)
        title('2D Deviation from Ground Truth, All Analyses','FontSize',22)
        % determine if all the outliers are removed
        outliersRemoved = input('Are all the cylinders in artifact regions removed? (y = 1, n = 0)');
    end
    % once all the outliers are removed plot the deviation
    plotDeviation(radiusFinalCombined,totDiffFinalCombined,binSize)
else
    plotDeviation(radiusAll,deviation2DAll,binSize)
end



% if you removed any cylinders, recalculate statistics, save new images and
% save the excel file
if totRemoved > 0
    cycleScans = 1;
    clear 'finalCircleDataX' 'finalCircleDataY' 'xGndTruthFinal' 'yGndTruthFinal'...
        'xCircle' 'yCircle' 'xCorrGndTruth' 'yCorrGndTruth' 'imageData'...
        'sliceShift' 'voxelLength'
    finalCircleDataX = xLocFinal{cycleScans};
    finalCircleDataY = yLocFinal{cycleScans};
    xGndTruthFinal = xGndFinal{cycleScans};
    yGndTruthFinal = yGndFinal{cycleScans};
    xCircle = xCircleFinal{cycleScans};
    yCircle = yCircleFinal{cycleScans};
    xCorrGndTruth = xCorrGndFinal{cycleScans};
    yCorrGndTruth = yCorrGndFinal{cycleScans};
    % non-upscaled image
    imageData = imageDataAll{cycleScans};
    sliceShift = sliceShiftExcel{cycleScans};
    voxelLength = voxelLengthExcel{cycleScans};
    
    
    % re-run procrustes on the analyses
    procylinderBefore = [finalCircleDataX',finalCircleDataY'];
    proGround = [xGndTruthFinal',yGndTruthFinal'];
    [d, procylinderAfter,transform] = procrustes(proGround,procylinderBefore,'scaling',false);
    proX = procylinderAfter(:,1);
    proY = procylinderAfter(:,2);
    finalCircleDataX = proX';
    finalCircleDataY = proY';
    
    % apply the opposite transformation on the ground truth
    proC = transform.c;
    proT = transform.T;
    % apply the reverse of procrustes --> Y = (Z - C)T^-1
    groundBeforeCorr = [xGndTruthFinal',yGndTruthFinal'];
    proCtemp = proC(1,:);
    xProC = proCtemp(1,1).*ones(length(xGndTruthFinal),1);
    yProC = proCtemp(1,2).*ones(length(xGndTruthFinal),1);
    proCsmall = [xProC,yProC];
    gndCorrected = (groundBeforeCorr - proCsmall)/(proT); % this does the inverse
    
    avgDiffGroundCorr = 0;
    avgDiffCylCorr = 0;
    for step = 1:length(gndCorrected(:,1))
        xCorrGndTruth(step) = gndCorrected(step,1);
        yCorrGndTruth(step) = gndCorrected(step,2);
        % after correcting the ground truth data instead
        avgDiffGroundCorr = avgDiffGroundCorr + sqrt((xCorrGndTruth(step)-xCircle(step)).^2 + ...
            (yCorrGndTruth(step)-yCircle(step)).^2);
        avgDiffCylCorr = avgDiffCylCorr + sqrt((finalCircleDataX(step)-xGndTruthFinal(step)).^2 + ...
            (finalCircleDataY(step)-yGndTruthFinal(step)).^2);
    end
    
    if sign(sliceShift) == 1
        posOrNeg = 'Pos';
    elseif sign(sliceShift) == 0
        posOrNeg = 'Iso';
    else
        posOrNeg = 'Neg';
    end
    nDataFinal = length(finalCircleDataX);
    deviation2D = zeros(1,nDataFinal);
    xDeviation = deviation2D;
    yDeviation = deviation2D;
    % calculate distance from center of the image
    radius = deviation2D;
    % check if dimensions are even or odd. If odd, just divide by 2 (center
    % position is an integer value), otherwise add 0.5
    if mod(length(imageData(1,:)),2) ~= 0
        xCenter = length(imageData(1,:))/2;
    else
        xCenter = length(imageData(1,:))/2 + 0.5;
    end
    if mod(length(imageData(1,:)),2) ~= 0
        yCenter = length(imageData(:,1))/2;
    else
        yCenter = length(imageData(:,1))/2 + 0.5;
    end
    for step = 1:length(xCircle)
        radius(step) = sqrt( (voxelWidth*(xCenter -  finalCircleDataX(step)))^2 + ...
            (voxelHeight*(yCenter -  finalCircleDataY(step)))^2 + sliceShift.^2);
        deviation2D(step) = sqrt( (voxelWidth*(finalCircleDataX(step) -  xGndTruthFinal(step)))^2 + ...
            (voxelHeight*(finalCircleDataY(step) -  yGndTruthFinal(step)))^2);
        xDeviation(step) =  sqrt( (voxelWidth*(finalCircleDataX(step) -  xGndTruthFinal(step)))^2);
        yDeviation(step) =  sqrt( (voxelHeight*(finalCircleDataY(step) -  yGndTruthFinal(step)))^2);
    end
    maxDevAll = max(deviation2D);
    [compareDeviationXY,avgDevAllXY,percentPassXY1,percentPassXY2,percentPassXY3,...
        avgDevRadius1XY,avgDevRadius2XY,avgDevRadius3XY,...
        maxDevRadius1XY,maxDevRadius2XY,maxDevRadius3XY]  =...
        compViewRay2D(finalCircleDataX,finalCircleDataY,xGndTruthFinal,yGndTruthFinal,...
        voxelHeight,voxelWidth,imageData,radiusThreshold1,radiusThreshold2,radiusThreshold3,...
        xCenter,yCenter,radiusSearch1,radiusSearch2,radiusSearch3,upscaleFactor,sliceShift,...
        xCircle,yCircle,xCorrGndTruth,yCorrGndTruth,plotData,plotOutsideRegion,MRorCT);
    
    if saveImages == 1;
        cd(saveImageDirectory)
        fig = gcf;
        fig.PaperPositionMode = 'auto';
        plotTitle = strcat(['compareViewRay',posOrNeg,num2str(abs(sliceShift)),...
            'SliceThick',num2str(round(voxelLength))]);
        print(plotTitle,'-dtiff','-r300')
        cd(currentDir)
        close(fig)
    end
    % plot deviation with respect to distance from center of the image
    if plotData == 1
        figure;
        plot(radius,deviation2D,'.')
        xlabel('Distance from Isosenter (mm)','FontSize',16)
        ylabel('2D Deviation (mm)','FontSize',16)
        title('2D Deviation from Gound Truth','FontSize',16)
    end
    if saveImages == 1;
        cd(saveImageDirectory)
        fig = gcf;
        fig.PaperPositionMode = 'auto';
        plotTitle = strcat(['IndDevPlotShift',posOrNeg,num2str(abs(sliceShift)),...
            'SliceThick',num2str(round(voxelLength))]);
        print(plotTitle,'-dtiff','-r300')
        cd(currentDir)
        close(fig)
    end
    
    % plot deviation in bins
    if plotData == 1
        plotDeviation(radius,deviation2D,binSize)
    end
    if saveImages == 1;
        cd(saveImageDirectory)
        fig = gcf;
        fig.PaperPositionMode = 'auto';
        plotTitle = strcat(['DevPlotShift',posOrNeg,num2str(abs(sliceShift)),...
            'SliceThick',num2str(round(voxelLength))]);
        print(plotTitle,'-dtiff','-r300')
        cd(currentDir)
        close(fig)
    end
    
    % store the data for saving the excel files
    percentPassXY1Excel{cycleScans} = percentPassXY1;
    percentPassXY2Excel{cycleScans} = percentPassXY2;
    percentPassXY3Excel{cycleScans} = percentPassXY3;
    avgDevRadius1XYExcel{cycleScans} = avgDevRadius1XY;
    avgDevRadius2XYExcel{cycleScans} = avgDevRadius2XY;
    avgDevRadius3XYExcel{cycleScans} = avgDevRadius3XY;
    avgDevAllXYExcel{cycleScans} = avgDevAllXY;
    maxDevRadius1XYExcel{cycleScans} = maxDevRadius1XY;
    maxDevRadius2XYExcel{cycleScans} = maxDevRadius2XY;
    maxDevRadius3XYExcel{cycleScans} = maxDevRadius3XY;
    maxDevAllExcel{cycleScans} = maxDevAll;
    radiusFinal{cycleScans} = radius;
    xLocFinal{cycleScans} = finalCircleDataX;
    yLocFinal{cycleScans} = finalCircleDataY;
    xGndFinal{cycleScans} = xGndTruthFinal;
    yGndFinal{cycleScans} = yGndTruthFinal;
    totDiffFinal{cycleScans} = deviation2D;
    % since we are only doing one scan
    percentPassAll = (countPass-totRemoved)/(nTotData--totRemoved)*100;
    percentPassAllExcel{cycleScans} = percentPassAll;
end

 
thisX = finalCircleDataX;
thisY = finalCircleDataY;
thisGndX = xGndTruthFinal;
thisGndY = yGndTruthFinal;
saveDeformationImg = 0;
devThisImg = zeros(1,length(thisX));
for step = 1:length(thisX)
    devThisImg(step) = sqrt((voxelWidth*(thisX(step) - thisGndX(step)))^ 2 + ...
        (voxelHeight*(thisY(step) - thisGndY(step)))^ 2);
end
plotDeformationField2D(thisX,thisY,thisGndX,thisGndY,devThisImg,upscaleImg,...
voxelWidth,voxelHeight,upscaleFactor,saveDeformationImg)

if saveExcelFiles == 1;
    saveExcel2D(corrThresholdExcel,radiusThreshold1Excel,radiusThreshold2Excel,...
        radiusThreshold3Excel,percentPassAllExcel,radiusSearch1Excel,radiusSearch2Excel,...
        radiusSearch3Excel,percentPassXY1Excel,percentPassXY2Excel,percentPassXY3Excel,avgDevRadius1XYExcel,avgDevRadius2XYExcel,...
        avgDevRadius3XYExcel,avgDevAllXYExcel,maxDevRadius1XYExcel,maxDevRadius2XYExcel,maxDevRadius3XYExcel,maxDevAllExcel,sliceShiftExcel,...
        voxelWidthExcel,voxelHeightExcel,voxelLengthExcel,fileNameExcel)
end