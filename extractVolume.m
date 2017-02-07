% Script to provide an example of how to extract images from DAT files
%
% Input:
% filename The name of the data file (.dat) you want to extract in the Raw Data folder
% fileDirectory The location of the file 
% maxSlice The number of image acquired in the dataset
% saveDicom (y = 1, n = 0) Whether or not you want to save the DICOM images
% 
% Output:
% averageImg The image resulting from the average of all the images
%
% John Ginn
% Created: 11/2/16
% Modified: 11/3/16

function [averageImg] = extractVolume(filename,fileDirectory,maxSlice,saveDicom)

% parameters
someDicomFileName = 'testing.dcm';
displayModulo = maxSlice; % Display input image every so often
outputDirectory = pwd; % The current directory
currentDirectory = pwd;
seriesNumber = 1;

% Create a dicom file for the function to overwrite
blankImg = zeros(40,40);
dicomwrite(blankImg,someDicomFileName)

% the series of images
[imgSeries] = ConvertSendExternalToDicom(filename, someDicomFileName, displayModulo, maxSlice, outputDirectory, seriesNumber,saveDicom,...
    fileDirectory,currentDirectory);

% average the images together
numOfImg = length(imgSeries(1,1,:));
averageImg = zeros(length(imgSeries(:,1,1)),length(imgSeries(1,:,1)));
for step = 1:numOfImg
    averageImg(:,:) = imgSeries(:,:,step); 
end
averageImg = averageImg./numOfImg;

close(gcf)
