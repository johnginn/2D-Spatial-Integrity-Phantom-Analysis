function varargout = RemoveCylinderGUI_single(varargin)
% REMOVECYLINDERGUI_SINGLE MATLAB code for RemoveCylinderGUI_single.fig
%      REMOVECYLINDERGUI_SINGLE, by itself, creates a new REMOVECYLINDERGUI_SINGLE or raises the existing
%      singleton*.
%
%      H = REMOVECYLINDERGUI_SINGLE returns the handle to a new REMOVECYLINDERGUI_SINGLE or the handle to
%      the existing singleton*.
%
%      REMOVECYLINDERGUI_SINGLE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in REMOVECYLINDERGUI_SINGLE.M with the given input arguments.
%
%      REMOVECYLINDERGUI_SINGLE('Property','Value',...) creates a new REMOVECYLINDERGUI_SINGLE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before RemoveCylinderGUI_single_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to RemoveCylinderGUI_single_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help RemoveCylinderGUI_single

% Last Modified by GUIDE v2.5 03-Feb-2017 17:06:39

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @RemoveCylinderGUI_single_OpeningFcn, ...
                   'gui_OutputFcn',  @RemoveCylinderGUI_single_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before RemoveCylinderGUI_single is made visible.
function RemoveCylinderGUI_single_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to RemoveCylinderGUI_single (see VARARGIN)

% Choose default command line output for RemoveCylinderGUI_single
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

if isappdata(hObject,'spheresToRemove')
    rmappdata(hObject,'spheresToRemove')
end

% obtain the data from the input cell data
xCircleFinal = varargin{1};
yCircleFinal = varargin{2};
xCorrGndFinal = varargin{3};
yCorrGndFinal = varargin{4};
imageDataAll = varargin{5};
sliceShiftExcel = varargin{6};
voxelLengthExcel = varargin{7};
totScans = varargin{8};
upscaleFactor = varargin{9};
MRorCT = varargin{10};
voxelWidth = varargin{11};
voxelHeight = varargin{12};
radiusSearch1 = varargin{13}; 
radiusSearch2 = varargin{14}; 
radiusSearch3 = varargin{15}; 
cycleScans = 1;

setappdata(hObject, 'countPlotFunctionCalls',0)
setappdata(hObject, 'countSpheresToRemove',0)
setappdata(hObject,'xCircleFinal',xCircleFinal);
setappdata(hObject,'yCircleFinal',yCircleFinal);
setappdata(hObject,'xCorrGndFinal',xCorrGndFinal);
setappdata(hObject,'yCorrGndFinal',yCorrGndFinal);
setappdata(hObject,'imageDataAll',imageDataAll);
setappdata(hObject,'sliceShiftExcel',sliceShiftExcel);
setappdata(hObject,'voxelLengthExcel',voxelLengthExcel);
setappdata(hObject,'cycleScans',cycleScans);
setappdata(hObject,'upscaleFactor',upscaleFactor);
setappdata(hObject,'MRorCT',MRorCT);
setappdata(hObject,'voxelWidth',voxelWidth);
setappdata(hObject,'voxelHeight',voxelHeight);
setappdata(hObject,'radiusSearch1',radiusSearch1);
setappdata(hObject,'radiusSearch2',radiusSearch2);
setappdata(hObject,'radiusSearch3',radiusSearch3);

disp('Close GUI when finished')

plotImage(handles)
% UIWAIT makes RemoveCylinderGUI_single wait for user response (see UIRESUME)
uiwait(handles.figure1);



function plotImage(handles)
% obtain data
countPlotFunctionCalls = getappdata(gcf,'countPlotFunctionCalls');
countPlotFunctionCalls = countPlotFunctionCalls + 1;
setappdata(gcf,'countPlotFunctionCalls',countPlotFunctionCalls);

xCircleFinal = getappdata(gcf,'xCircleFinal');
yCircleFinal = getappdata(gcf,'yCircleFinal');
xCorrGndFinal = getappdata(gcf,'xCorrGndFinal');
yCorrGndFinal = getappdata(gcf,'yCorrGndFinal');
imageDataAll = getappdata(gcf,'imageDataAll');
sliceShiftExcel = getappdata(gcf,'sliceShiftExcel');
voxelLengthExcel = getappdata(gcf,'voxelLengthExcel');
cycleScans = getappdata(gcf,'cycleScans');
upscaleFactor = getappdata(gcf,'upscaleFactor');
MRorCT = getappdata(gcf,'MRorCT');
voxelWidth = getappdata(gcf,'voxelWidth');
voxelHeight = getappdata(gcf,'voxelHeight');
radiusSearch1 = getappdata(gcf,'radiusSearch1');
radiusSearch2 = getappdata(gcf,'radiusSearch2');
radiusSearch3 = getappdata(gcf,'radiusSearch3');



% extract the current scan data
xCircle = xCircleFinal{cycleScans};
yCircle = yCircleFinal{cycleScans};
xCorrGndTruth = xCorrGndFinal{cycleScans};
yCorrGndTruth = yCorrGndFinal{cycleScans};
% non-upscaled image
imageData = imageDataAll{cycleScans};
sliceShift = sliceShiftExcel{cycleScans};
voxelLength = voxelLengthExcel{cycleScans};

if sign(sliceShift) == 1
    posOrNeg = 'Pos';
elseif sign(sliceShift) == 0
    posOrNeg = 'Iso';
else
    posOrNeg = 'Neg';
end
centerRow = round(median(1:length(imageData(:,1))));
centerCol = round(median(1:length(imageData(1,:))));

cylinderImg = zeros(upscaleFactor*length(imageData(:,1)),upscaleFactor*length(imageData(1,:)));
groundImg = zeros(upscaleFactor*length(imageData(:,1)),upscaleFactor*length(imageData(1,:)));
for step = 1:length(xCircle)
    data = round([yCircle(step) xCircle(step)].*upscaleFactor);
    xMarkerLoc = (data(1)-upscaleFactor-1):1:(data(1)+upscaleFactor+1); % make + symbol on image
    yMarkerLoc = (data(2)-upscaleFactor-1):1:(data(2)+upscaleFactor+1);% make + symbol on image
    cylinderImg(xMarkerLoc,yMarkerLoc) = makeShape(length(xMarkerLoc),'square','thick'); % make a square
    
    data = round([yCorrGndTruth(step) xCorrGndTruth(step)].*upscaleFactor);
    xMarkerLoc = (data(1)-upscaleFactor):1:(data(1)+upscaleFactor); % make + symbol on image
    yMarkerLoc = (data(2)-upscaleFactor):1:(data(2)+upscaleFactor);% make + symbol on image
    groundImg(xMarkerLoc,data(2)) = 1; % - portion of + marker for ground truth of cylinders
    groundImg(data(1),yMarkerLoc) = 1; % - portion of + marker for ground truth of cylinders
    groundImg(xMarkerLoc,data(2)-1) = 1; % - portion of + marker for ground truth of cylinders
    groundImg(data(1)-1,yMarkerLoc) = 1; % - portion of + marker for ground truth of cylinders
    groundImg(xMarkerLoc,data(2)+1) = 1; % - portion of + marker for ground truth of cylinders
    groundImg(data(1)+1,yMarkerLoc) = 1; % - portion of + marker for ground truth of cylinders
end
% grab the handle of the first plot
axes(handles.Plot1)

if (~strcmp(MRorCT,'ct'))
    % upscale the image for plotting
    % original information
    origX = 1:length(imageData(1,:));
    origY = 1:length(imageData(:,1));
    [origMeshX,origMeshY] = meshgrid(origX,origY);
    % new data
    volUpscaleX = 1/upscaleFactor:1/upscaleFactor:length(imageData(1,:));
    volUpscaleY = 1/upscaleFactor:1/upscaleFactor:length(imageData(:,1));
    [upMeshX, upMeshY] = meshgrid(volUpscaleX,volUpscaleY);
    % the upscaled data
    SigDataUpscale = interp2(origMeshX,origMeshY,imageData,...
        upMeshX,upMeshY,'cubic');
    
    upVoxWidth = voxelWidth/upscaleFactor;
    upVoxHeight = voxelHeight/upscaleFactor;
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
    SigDataUpscale = imageData;
    upCenterRow = centerRow;
    upCenterCol = centerCol;
    upVoxWidth = voxelWidth;
    upVoxHeight = voxelHeight;
end



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
green = cat(3, zeros(size(cylinderImg)),ones(size(cylinderImg)), zeros(size(cylinderImg)));
% bright = 0.9;
% white version
% green = bright.*cat(3, ones(size(cylinderFitPtsPass)),ones(size(cylinderFitPtsPass)), ones(size(cylinderFitPtsPass)));
hold on
hGreen = imshow(green);
hold off
set(hGreen, 'AlphaData', cylinderImg) % make color sheet only show markers

% ground truth marker plotting
% cat(3,r,g,b)
white = cat(3, ones(size(groundImg)),ones(size(groundImg)),ones(size(groundImg)));
hold on
hWhite = imshow(white);
hold off
set(hWhite, 'AlphaData', groundImg) % make color sheet only show markers
% title('Spatial Integrity Phantom Center Slice','FontSize',20)

if countPlotFunctionCalls == 1
   % this is the first time this function is called after creating GUI object
   plotRemove = zeros(length(groundImg(:,1)),length(groundImg(1,:)));
   setappdata(gcf,'plotRemove',plotRemove)
end
plotRemove = getappdata(gcf,'plotRemove');
red = cat(3, ones(size(plotRemove)),zeros(size(plotRemove)),zeros(size(plotRemove)));
hold on
hWhite = imshow(red);
hold off
set(hWhite, 'AlphaData', plotRemove) % make color sheet only show markers

xlabel('y-position (mm)','FontSize',20)
ylabel('z-position (mm)','FontSize',20)
title(strcat(['Cylinder Deviation Test, ',num2str(sliceShift),' (mm) shift']),'FontSize',20)



% --- Outputs from this function are returned to the command line.
function varargout = RemoveCylinderGUI_single_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
% varargout{1} = handles.output;
varargout{1} = getappdata(gcf,'spheresToRemove');

delete(handles.figure1)



% --- Executes on button press in markButton.
function markButton_Callback(hObject, eventdata, handles)
% hObject    handle to markButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

upscaleFactor = getappdata(gcf,'upscaleFactor');
countSpheresToRemove = getappdata(gcf,'countSpheresToRemove');
plotRemove = getappdata(gcf,'plotRemove');
scanNumber = getappdata(gcf,'cycleScans');
% if this isn't the first sphere to be removed obtain the others in the
% list
if countSpheresToRemove ~= 0
    spheresToRemove = getappdata(gcf,'spheresToRemove');
end
countSpheresToRemove = countSpheresToRemove + 1;
setappdata(gcf,'countSpheresToRemove',countSpheresToRemove)
% obtain the data locations
[xRemove, yRemove] = ginput(1);
spheresToRemove(countSpheresToRemove,1) = xRemove;
spheresToRemove(countSpheresToRemove,2) = yRemove;
spheresToRemove(countSpheresToRemove,3) = scanNumber;
xDim = length(plotRemove(1,:));
yDim = length(plotRemove(:,1));
xPos = round(xRemove);
yPos = round(yRemove);
markerSize = 3; % size on each side of marker
if (yPos > 1)&&(yPos < yDim)
    yLoc = (yPos-markerSize):1:(yPos+markerSize);
end
if (yPos <= 1)
    yLoc = 1:1:(2);
end
if (yPos >= yDim)
    yLoc = (yDim-1):1:yDim;
end
% x-location
if (xPos > 1)&&(xPos < xDim)
    xLoc = (xPos-markerSize):1:(xPos+markerSize);
end
if (xPos <= 1)
    xLoc = 1:1:(2);
end
if (xPos >= xDim)
    xLoc = (xDim-1):1:xDim;
end
plotRemove(yLoc,xLoc) = 1;
setappdata(gcf,'spheresToRemove',spheresToRemove);
setappdata(gcf,'plotRemove',plotRemove)
plotImage(handles);


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
if isequal(get(hObject, 'waitstatus'),'waiting')
    uiresume(hObject);
else
    delete(hObject);
end
