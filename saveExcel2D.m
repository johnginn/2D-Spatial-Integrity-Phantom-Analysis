% This function saves the different variables obtained from the analysis to
% excel files. WARNING: Be sure to transfer over data from a previous
% analysis prior to running this function, otherwise the data will be
% overwritten
%
% Input:
% corrThreshold The correlation coefficient threshold to determine if the cylinder was found
% radiusThreshold1 The deviation threshold defined for cylinders within radiusSearch1
% radiusThreshold2 The deviation threshold defined for cylinders within radiusSearch2
% radiusThreshold3 The deviation threshold defined for cylinders within radiusSearch3
% percentPassAll The percent of all cylinders passing the correlation threshold
% radiusSearch1 The radius from isocenter defining the first search region
% radiusSearch2 The radius from isocenter defining the second search region
% radiusSearch3 The radius from isocenter defining the third search region
% percentPassXY1 The percentage of cylinders passing radiusThreshold1 within radiusSearch1
% percentPassXY2 The percentage of cylinders passing radiusThreshold2 within radiusSearch2
% percentPassXY3 The percentage of cylinders passing radiusThreshold3 within radiusSearch3
% avgDevRadius1XY The average deviation for all cylinders within radiusSearch1
% avgDevRadius2XY The average deviation for all cylinders within radiusSearch2
% avgDevRadius3XY The average deviation for all cylinders within radiusSearch3
% avgDevAllXY  The average deviation for all cylinders
% maxDevRadius1XY The maximum deviation for all cylinders within radiusSearch1
% maxDevRadius2XY The maximum deviation for all cylinders within radiusSearch2
% maxDevRadius3XY The maximum deviation for all cylinders within radiusSearch3
% maxDevAll The maximum deviation for all cylinders
% sliceShift (mm) The shift away from isocenter in the lateral direction
% voxelWidith (mm) The width of the voxel x-dimension (medial-lateral)
% voxelHeight (mm) The width of the voxel y-dimension (anterior-posterior)
% voxelLength (mm) The width of the voxel z-dimension (inferior-superior)
% fileName The name of the file containing the real-time imaging data
%
% Output:
%
% John Ginn
% Created: 11/8/16
% Modified: 11/8/16

function [] = saveExcel2D(corrThreshold,radiusThreshold1,radiusThreshold2,...
        radiusThreshold3,percentPassAll,radiusSearch1,radiusSearch2,...
        radiusSearch3,percentPassXY1,percentPassXY2,percentPassXY3,avgDevRadius1XY,avgDevRadius2XY,...
        avgDevRadius3XY,avgDevAllXY,maxDevRadius1XY,maxDevRadius2XY,maxDevRadius3XY,maxDevAll,sliceShift,...
        voxelWidth,voxelHeight,voxelLength,fileName)

% table containing the data
saveData = cell(15,(length(corrThreshold) + 1));
for stepData = 1:(length(corrThreshold) + 1)
    % just store the titles for reading the data
    if stepData == 1;
        saveData{1,stepData} = 'File Name:';
        saveData{2,stepData} = 'Shift from isocenter';
        saveData{3,stepData} = 'Resolution';
        saveData{4,stepData} = ...
            strcat(['Percent of cylinders passing correlation threshold of ',num2str(corrThreshold{stepData})]);
        saveData{5,stepData} = strcat(['Cylinders w/in ',num2str(radiusSearch1{stepData}),' (mm) of image center,'...
            ' less than ',num2str(radiusThreshold1{stepData}),' (mm) deviation']);
        saveData{6,stepData} = strcat(['Cylinders w/in ',num2str(radiusSearch2{stepData}),' (mm) of image center,'...
            ' less than ',num2str(radiusThreshold2{stepData}),' (mm) deviation']);
        saveData{7,stepData} = strcat(['Cylinders w/in ',num2str(radiusSearch3{stepData}),' (mm) of image center,'...
            ' less than ',num2str(radiusThreshold3{stepData}),' (mm) deviation']);
        saveData{8,stepData} = 'Avg 2D deviation for all cylinders';
        saveData{9,stepData} = strcat(['Cylinders w/in ',num2str(radiusSearch1{stepData}),' (mm) of image center,'...
            ' avg 2D deviation']);
        saveData{10,stepData} = strcat(['Cylinders w/in ',num2str(radiusSearch2{stepData}),' (mm) of image center,'...
            ' avg 2D deviation']);
        saveData{11,stepData} = strcat(['Cylinders w/in ',num2str(radiusSearch3{stepData}),' (mm) of image center,'...
            ' avg 2D deviation']);
        saveData{12,stepData} = 'Max 2D deviation for all cylinders';
        saveData{13,stepData} = strcat(['Max 2D deviation for cylinders w/in ',num2str(radiusSearch1{stepData}),' (mm) of image center']);
        saveData{14,stepData} = strcat(['Max 2D deviation for cylinders w/in ',num2str(radiusSearch2{stepData}),' (mm) of image center']);
        saveData{15,stepData} = strcat(['Max 2D deviation for cylinders w/in ',num2str(radiusSearch3{stepData}),' (mm) of image center']);
    else
        saveData{1,stepData} = fileName{stepData - 1};
        saveData{2,stepData} = strcat([num2str(sliceShift{stepData - 1}),' (mm)']);
        saveData{3,stepData} = strcat([num2str(voxelWidth{stepData - 1}),'x',...
            num2str(voxelHeight{stepData - 1}),'x',num2str(voxelLength{stepData - 1}),' (mm)']);...
            saveData{4,stepData} = strcat([num2str(percentPassAll{stepData - 1}),' (%)']);
        saveData{5,stepData} = strcat([num2str(percentPassXY1{stepData - 1}),' (%)']);
        saveData{6,stepData} = strcat([num2str(percentPassXY2{stepData - 1}),' (%)']);
        saveData{7,stepData} = strcat([num2str(percentPassXY3{stepData - 1}),' (%)']);
        saveData{8,stepData} = strcat([num2str(avgDevAllXY{stepData - 1}),' (mm)']);
        saveData{9,stepData} = strcat([num2str(avgDevRadius1XY{stepData - 1}),' (mm)']);
        saveData{10,stepData} = strcat([num2str(avgDevRadius2XY{stepData - 1}),' (mm)']);
        saveData{11,stepData} = strcat([num2str(avgDevRadius3XY{stepData - 1}),' (mm)']);
        saveData{12,stepData} = strcat([num2str(maxDevAll{stepData - 1}),' (mm)']);
        saveData{13,stepData} = strcat([num2str(maxDevRadius1XY{stepData - 1}),' (mm)']);
        saveData{14,stepData} = strcat([num2str(maxDevRadius2XY{stepData - 1}),' (mm)']);
        saveData{15,stepData} = strcat([num2str(maxDevRadius3XY{stepData - 1}),' (mm)']);
    end
end
filename = 'AnalysisResults.xlsx';
xlswrite(filename,saveData);

end