% Function to read the filenames for the analysis
%
% Input:
%
% Output:
% filenameSorted The filenames sorted according to the time they were acquired
% sortedTime The time the image was acquired
% fileDirectory The location of the file
%
% John Ginn
% Created: 11/9/16
% Modified: 11/9/16

function [filenameSorted,sortedTime,fileDirectory] = readFilenames()
currentDir = pwd;
disp('Please select folder containing the files to be converted to DICOM')
fileDirectory = uigetdir;
cd(fileDirectory)
filesStructure = dir;
cd(currentDir);

count = 0; % then number of data files
% sort through files and use only .dat files
for step = 1:(length(filesStructure) - 2)
    % +2 to ignore the first two lines '.' and '..' respectively
    currentFile = filesStructure(step + 2).name;
    if strcmp('.dat',currentFile((length(currentFile)-3):length(currentFile)))
        count = count + 1;
        fileNames{count} = currentFile;
    end
    
end
timeArray = zeros(length(fileNames),1);
for step = 1:length(fileNames);
    currentString = fileNames{step};
    % determine the time of each scan based on the file name
    currStringLeng = length(currentString);
    hour = str2num(currentString((currStringLeng-19):(currStringLeng-18)));
    min = str2num(currentString((currStringLeng-16):(currStringLeng-15)));
    sec = str2num(currentString((currStringLeng-13):(currStringLeng-12)));
    timeArray(step) = hour*60*60 + min*60 + sec;
end

% sort the times of the scans
[sortedTime, sortedI] = sort(timeArray);

% sort the filenames
filenameSorted = cell(length(fileNames),1);
for step = 1:length(fileNames)
    filenameSorted{step} = fileNames{sortedI(step)};
end




end