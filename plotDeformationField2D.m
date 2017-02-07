% Function to plot the deformation field map for a slice in the phantom.
% NOTE: the indices specified below MUST correspond to the locations in the
% image. Additinoally, the ground truth locations must not be corrected
% by a transform
%
% Input:
% xSph (index) The sphere found x-location
% ySph (index) The sphere found y-location
% xGnd (index) The ground truth x-location
% yGnd (index) The ground truth y-location
% deviation (mm) The deviation between the sphere and ground truth locations
% sliceImg The image of the slice for these spheres 
% voxelWidth (mm) The width of the voxels
% voxelHeight (mm) The height of the voxels
% scaleFactor The factor the image is upsampled by
% saveDeformationImg (y=1,n=0) Whether or not to save the deformation field image
%
% Output:
%
% John Ginn
% Created: 1/25/17
% Modified: 1/30/17

function [] = plotDeformationField2D(xSph,ySph,xGnd,yGnd,deviation,sliceImg,...
    voxelWidth,voxelHeight,scaleFactor,saveDeformationImg)
scaleMin = 0;
scaleMax = 3;
extraDist = 3; % 3 mm extra plotting


xImgDim = length(sliceImg(1,:));
yImgDim = length(sliceImg(:,1));
% extra distance for plotting
% find bounds of region for deformation map
minX = round((min(xSph) - extraDist/voxelWidth*scaleFactor)*scaleFactor);
maxX = round((max(xSph) + extraDist/voxelWidth*scaleFactor)*scaleFactor);
minY = round((min(ySph) - extraDist/voxelHeight*scaleFactor)*scaleFactor);
maxY = round((max(ySph) + extraDist/voxelHeight*scaleFactor)*scaleFactor);
if minX < 1
   minX = 1; 
end
if maxX > xImgDim
   maxX = xImgDim; 
end
if minY < 1
    minY = 1;
end
if maxY > yImgDim
   maxY = yImgDim; 
end

xArray = minX:maxX;
yArray = minY:maxY;
[meshX,meshY] = meshgrid(xArray,yArray);


minGndX = (min(xGnd)).*scaleFactor;
maxGndX = (max(xGnd)).*scaleFactor;
minGndY = (min(xGnd)).*scaleFactor;
maxGndY = (max(ySph)).*scaleFactor;

xGndUnique = unique(xGnd.*scaleFactor);
yGndUnique = unique(yGnd.*scaleFactor);

[meshGndX,meshGndY] = meshgrid(xGndUnique,yGndUnique);
devMesh = zeros(length(meshGndX(:,1)),length(meshGndX(1,:)));
for step = 1:length(xSph)
    % find the location of this value in the grid
    thisX = find(xGnd(step).*scaleFactor == xGndUnique);
    thisY = find(yGnd(step).*scaleFactor == yGndUnique);
    devMesh(thisY,thisX) = deviation(step);
end
X = meshGndX;
Y = meshGndY;
V = devMesh;
Xq = meshX;
Yq= meshY;
interpVal = interp2(X,Y,V,Xq,Yq,'linear',0);
% remove any deformation less than 0
for stepX = 1:length(interpVal(1,:))
    for stepY = 1:length(interpVal(:,1))
        if interpVal(stepY,stepX) < scaleMin;
            interpVal(stepY,stepX) = 0;
        end
    end
end
% create the interpolated image. interpVal is not necessarily the same dimension as
% the image
interpDefImg = zeros(length(sliceImg(:,1)),length(sliceImg(1,:)));
interpDefImg(yArray,xArray) = interpVal;
% create the mask for this image
interpMask = zeros(length(sliceImg(:,1)),length(sliceImg(1,:)));
for stepX = 1:length(interpDefImg(1,:))
    for stepY = 1:length(interpDefImg(:,1))
        if interpDefImg(stepY,stepX) > scaleMin;
            interpMask(stepY,stepX) = 0.5;
        end
    end
end
% normalize the image scale of the phantom to the deformation map
maxDef = max(max(interpDefImg));
minDef = 0;


normImg = sliceImg;
minVal = min(min(normImg));
% make sure the smallest value is at least zero
if minVal < 0
    normImg(:) = normImg(:) + abs(minVal);
end
% normalize the image scale to the deformation scale
maxImg = max(max(normImg));
normFact = maxDef/maxImg;
normImg = normImg.*normFact; 

% shift the minimum value of the image to zero
minVal = min(min(normImg)); 
normImg = normImg - minVal;

% remove any NaN values from the image and set to zero (around the
% boundaries)
countNan = 0;
for stepX = 1:length(normImg(1,:))
    for stepY = 1:length(normImg(:,1))
        TF = isnan(normImg(stepY,stepX));
        if TF == 1
           normImg(stepY,stepX) = 0; 
           countNan = countNan + 1;
        end
    end
end


figure;
h = imagesc(normImg);
colormap('gray')
% Now make an RGB image that matches display from IMAGESC:
C = colormap;  % Get the figure's colormap.
L = size(C,1);
% Scale the matrix to the range of the map.
Gs = round(interp1(linspace(min(normImg(:)),max(normImg(:)),L),1:L,normImg));
colorPhantomImg = reshape(C(Gs,:),[size(Gs) 3]); % Make RGB image from scaled.
close gcf

% here is the plot that actually gets saved
figure;
imshow(colorPhantomImg,[]);
hold on
h = imshow(interpDefImg,[scaleMin scaleMax]);
hold off
colormap('jet')
barObj = colorbar('eastoutside');
% alpha 0.3
set(h, 'AlphaData', interpMask);
% change the colorbar ticks to include mm
currTicks = barObj.Ticks;
newTicks = cell(1,length(currTicks));
for step = 1:length(currTicks)
    newTicks{step} = strcat([num2str(currTicks(step)),' (mm)']);
end
newTicks{step} = strcat('>',[num2str(currTicks(step)),' (mm)']);
newTicks{1} = strcat([num2str(currTicks(1)),' (mm)']);
colorbar('Ticks',currTicks,...
         'TickLabels',newTicks)
% % This works with the exception of the colorbar
% figure;
% imshow(colorInterpImg);
% hold on
% h = imshow(normImg,[]);
% hold off
% % alpha 0.3
% set(h, 'AlphaData', interpMask);
% set(h,'colorbar',barObjTemp)
% barObj = colorbar('eastoutside');

% custom axes to show distances (the locations of the axes)
numOfTicks = 10;
xAxisSpacing = floor(length(sliceImg(1,:))/numOfTicks);
yAxisSpacing = floor(length(sliceImg(:,1))/numOfTicks);

% xAxisLocation = linspace(1,length(SigDataUpscale(1,:)),numOfTicks);
% yAxisLocation = linspace(1,length(SigDataUpscale(:,1)),numOfTicks);
xAxisLocation = 1:xAxisSpacing:(xAxisSpacing*numOfTicks+1);
yAxisLocation = 1:yAxisSpacing:(yAxisSpacing*numOfTicks+1);
% the values on the axes
xAxisValue = voxelWidth/scaleFactor.*(xAxisLocation - round(median(xAxisLocation)));
yAxisValue = voxelHeight/scaleFactor.*(yAxisLocation - round(median(yAxisLocation)));
xAxisLocation(length(xAxisLocation)) = length(sliceImg(1,:)); % special case for plotting for paper
yAxisLocation(length(yAxisLocation)) = length(sliceImg(:,1)); % special case for plotting for paper
axis on
axisHandle = gca;
set(gca,'XTickMode','manual')
set(gca,'YTickMode','manual')
set(gca,'XTick',xAxisLocation)
set(gca,'YTick',yAxisLocation)
set(gca,'XTickLabel',xAxisValue)
set(gca,'YTickLabel',yAxisValue)
set(gca,'FontSize',16)
xlabel('x-position (mm)','FontSize',20)
ylabel('z-position (mm)','FontSize',20)

if saveDeformationImg
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    print('DeformationFieldMap','-dtiff','-r300')
end
end