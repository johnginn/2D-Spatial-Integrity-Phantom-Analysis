% This function extracts file names of images from a folder containing 
% the dicom images exported from MIM of the spatial integrity phantom. The
% images are sorted based on their slice position.
%
% Input:
% folderDirectory The directory of the folder containing the images
% currentDirectory The directory where the script running the analysis is located
% linesToSkip (optional, do not use if unsure of lines to skip) Avoids extra lines when parsing the dicom file in the folder
% 
% Output:
% fileNames A cell array containing the strings of the file names
% fileData The image data from the dicom file
% fileMapThe color map of the dicom file
% fileInfo The dicom header information from the files
% 
% John Ginn
% created: 6/24/16
% modified: 8/16/16

function [fileNames, fileData, fileMap, fileInfo] =...
    loadImages(folderDirectory,currentDirectory,linesToSkip)
% change directory to location of folder conaining images
cd(folderDirectory)
% pull out file names
imFolderNames = dir;
% if linesToSkip does not exist
if (exist('linesToSkip') == 0)
    disp('skip lines until you see an image filename')
    % determine how many lines to skip because of odd name storing
    for step = 1:length(imFolderNames)
        disp(strcat('current line: ',imFolderNames(step).name))
        decision = input('skip this line? (y = 1/n = 0)');
        linesToSkip = step;
        if decision == 0;
            break
        end
    end
end
firstLine = linesToSkip + 1;

count = 1; % for storing the fileNames

for step = (firstLine - 1):length(imFolderNames) % weird issue with matlab -1
    fileNames{count,1} = imFolderNames(step).name;
    [fileData{count,1} fileMap{count,1}]= dicomread(fileNames{count,1});
    % convert to type double
    fileData{count} = im2double(fileData{count});
    fileInfo{count,1} = dicominfo(fileNames{count,1});
    count = count + 1;
end




% change back to location of code
cd(currentDirectory)

end