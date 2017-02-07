function [] = DicomToSendExternalSimpleCine( dicomPath, outputDirectory, outputFilename, numSlicesPerAcq )
%%
% Read dicom, write send external
% dicomPath where are the dicom files
% planningDirectory - output location for send external

display = 1;
logse = 1;

it=dir([ dicomPath '*']);    % look for existing cases

nit = length(it);
if nit == 0
    disp( 'Empty directory' );
    return;
end

% Open send external filename
filenameSendExternal = sprintf( '%s\\%s.dat', outputDirectory, outputFilename );
fidSendExternal = fopen( filenameSendExternal, 'w' );

%sortedPosVec = SortDicomDirFiles(dicomPath);

iAcqSliceIndex = 1;
%for ifile = 1:1:nit-2
for ifile = 3:1:nit
    image = dicomread( [dicomPath it(ifile).name] );
    %image = dicomread( [dicomPath it(sortedPosVec(ifile,1)).name] );
    image = permute( image, [ 2 1 ] );
    info = dicominfo( [dicomPath it(ifile).name] );
    %info = dicominfo( [dicomPath it(sortedPosVec(ifile,1)).name] );

    disp( sprintf( '%03d %02d %02d %02d %02d %02d %02d %10.5f %10.5f %10.5f ', ...
        ifile,...
        info.ImageOrientationPatient(4),...
        info.ImageOrientationPatient(5),...
        info.ImageOrientationPatient(6),...        
        info.ImageOrientationPatient(1),...
        info.ImageOrientationPatient(2),...
        info.ImageOrientationPatient(3),...        
        info.ImagePositionPatient(1),...        
        info.ImagePositionPatient(2),...        
        info.ImagePositionPatient(3) ) );

    header = BuildSendExternalHeader( info, iAcqSliceIndex );

    nsize = size(header);

    headerSize = nsize(2);
    counths = fwrite( fidSendExternal, headerSize, 'int32' );
    imageSize = double(info.Rows)*double(info.Columns)*2;
    countis = fwrite( fidSendExternal, imageSize, 'int32' );    
    counth = fwrite( fidSendExternal, header, 'char' );
    counti = fwrite( fidSendExternal, image, 'int16' );
    
    iAcqSliceIndex = iAcqSliceIndex + 1;
    if iAcqSliceIndex > numSlicesPerAcq
        iAcqSliceIndex = 1;
    end
end

disp( sprintf( 'Writing end record' ) );
zero = int32(0);
counths = fwrite( fidSendExternal, zero, 'int32' );
countis = fwrite( fidSendExternal, zero, 'int32' );
fclose( fidSendExternal );

end


function [sortedPosVec] = SortDicomDirFiles( dicomPath )
%%
it=dir([ dicomPath '*']);    % look for existing cases

nit = length(it);
if nit == 0
    disp( 'Empty directory!!!' );
    return;
end

% Use DICOM image info from first two images
% to see which dimension is changing for the sort.
% Image1 = file index 3 because . and .. are in the
% file list at index 1 and 2
info1 = dicominfo( [dicomPath it(3).name] );
info2 = dicominfo( [dicomPath it(4).name] );

disp( sprintf( 'Image 1 info %02d %02d %02d %02d %02d %02d %10.5f %10.5f %10.5f ', ...
    info1.ImageOrientationPatient(4),...
    info1.ImageOrientationPatient(5),...
    info1.ImageOrientationPatient(6),...        
    info1.ImageOrientationPatient(1),...
    info1.ImageOrientationPatient(2),...
    info1.ImageOrientationPatient(3),...        
    info1.ImagePositionPatient(1),...        
    info1.ImagePositionPatient(2),...        
    info1.ImagePositionPatient(3) ) );

disp( sprintf( 'Image 2 info %02d %02d %02d %02d %02d %02d %10.5f %10.5f %10.5f ', ...
    info2.ImageOrientationPatient(4),...
    info2.ImageOrientationPatient(5),...
    info2.ImageOrientationPatient(6),...        
    info2.ImageOrientationPatient(1),...
    info2.ImageOrientationPatient(2),...
    info2.ImageOrientationPatient(3),...        
    info2.ImagePositionPatient(1),...        
    info2.ImagePositionPatient(2),...        
    info2.ImagePositionPatient(3) ) );

% Determine which position dimension is changing for the sort
chgDim = 0;
if info1.ImagePositionPatient(1) ~= info2.ImagePositionPatient(1)
    chgDim = 1;
else
    if info1.ImagePositionPatient(2) ~= info2.ImagePositionPatient(2)
        chgDim = 2;
    else
        if info1.ImagePositionPatient(3) ~= info2.ImagePositionPatient(3)
            chgDim = 3;
        else
            disp( 'Could not determine sort order!!!' );
            return;
        end
    end
end

% Get posVecs for all files to prepare for the sor
posVec = zeros([nit-2,2]);
for ifile = 3:1:nit
    info = dicominfo( [dicomPath it(ifile).name] );

    disp( sprintf( '%03d %02d %02d %02d %02d %02d %02d %10.5f %10.5f %10.5f ', ...
        ifile,...
        info.ImageOrientationPatient(4),...
        info.ImageOrientationPatient(5),...
        info.ImageOrientationPatient(6),...        
        info.ImageOrientationPatient(1),...
        info.ImageOrientationPatient(2),...
        info.ImageOrientationPatient(3),...        
        info.ImagePositionPatient(1),...        
        info.ImagePositionPatient(2),...        
        info.ImagePositionPatient(3) ) );
    
    posVec(ifile-2,1) = ifile;
    posVec(ifile-2,2) = info.ImagePositionPatient(chgDim);
end

sortedPosVec = sortrows(posVec, 2);
end


function [posvec] = CreatePosVec( info )
%%
%
posvec = [ 0 0 0 ];

% Compute Posvec 0
nd2 = 0; 
spacing = 0;
if info.ImageOrientationPatient(4) == 1
    % colvec is x and axis is increasing
    nd2 = double(info.Rows)/2;
    spacing = info.PixelSpacing(1);
elseif info.ImageOrientationPatient(4) == -1 
    % colvec is x and axis is decreasing
    nd2 = double(info.Rows)/2;
    spacing = -info.PixelSpacing(1);
elseif info.ImageOrientationPatient(1) == 1
    % rowvec is y and axis is increasing
    nd2 = double(info.Columns)/2;
    spacing = info.PixelSpacing(2);
elseif info.ImageOrientationPatient(1) == -1 
    % colvec is y and axis is increasing
    nd2 = double(info.Columns)/2;
    spacing = -info.PixelSpacing(2);
end

posvec(1) = info.ImagePositionPatient(1) + nd2 * spacing;

% Posvec 1
nd2 = 0;
spacing = 0;
if info.ImageOrientationPatient(5) == 1
    % colvec is y and axis is increasing
    nd2 = double(info.Rows)/2;
    spacing = info.PixelSpacing(1);
elseif info.ImageOrientationPatient(5) == -1 
    % colvec is y and axis is decreasing
    nd2 = double(info.Rows)/2;
    spacing = -info.PixelSpacing(1);
elseif info.ImageOrientationPatient(2) == 1
    % rowvec is y and axis is increasing
    nd2 = double(info.Columns)/2;
    spacing = info.PixelSpacing(2);
elseif info.ImageOrientationPatient(2) == -1 
    % rowvec is y and axis is decreasing
    nd2 = double(info.Columns)/2;
    spacing = -info.PixelSpacing(2);
end

posvec(2) = info.ImagePositionPatient(2) + nd2 * spacing;
%posvec(2) = 0.0;

% Posvec 2
nd2 = 0;
spacing = 0;
if info.ImageOrientationPatient(6) == 1
    % colvec is z and axis is increasing
    nd2 = double(info.Columns)/2;
    spacing = info.PixelSpacing(1);
elseif info.ImageOrientationPatient(6) == -1 
    % colvec is z and axis is decreasing
    nd2 = double(info.Columns)/2;
    spacing = -info.PixelSpacing(1);
elseif info.ImageOrientationPatient(3) == 1
    % rowvec is z and axis is decreasing
    nd2 = double(info.Rows)/2;
    spacing = info.PixelSpacing(2);
elseif info.ImageOrientationPatient(3) == -1 
    % rowvec is z and axis is decreasing
    nd2 = double(info.Rows)/2;
    spacing = -info.PixelSpacing(2);
end

posvec(3) = info.ImagePositionPatient(3) + nd2 * spacing;

end
function [header] = BuildSendExternalHeader( info, sliceIndex )
%%
%
anatomicalSliceNumber = sliceIndex-1;
sliceNumber = sliceIndex-1;
%
% header1
s1 = 'ViewRay (c) Demo, Simulation Treatment data';
s2 = 'CONTROL.AnatomicalSliceNo = ';%0
s3 = 'CONTROL.ChronSliceNo = ';%1
s4 = 'DICOM.SliceNo = ';%0
s5 = 'DICOM.NoOfCols = ';%320
s6 = 'DICOM.NoOfRows = ';%270

%header2
s7 = 'DICOM.SliceThickness = ';%3.500000
s8 = 'DICOM.SliceLocation = ';%-199.752430
s9 = 'DICOM.PixelSpacing.0 = ';%1.296875
s10 = 'DICOM.PixelSpacing.1 = ';%1.296875

%header3
s11 = 'DICOM.ColVec.0 = ';%0.0
s12 = 'DICOM.ColVec.1 = ';%1.0
s13 = 'DICOM.ColVec.2 = ';%0.0
s14 = 'DICOM.RowVec.0 = ';%1.0
s15 = 'DICOM.RowVec.1 = ';%0.0
s16 = 'DICOM.RowVec.2 = ';%0.0
s17 = 'DICOM.PosVec.0 = ';%0.000000
s18 = 'DICOM.PosVec.1 = ';%13.320000
s19 = 'DICOM.PosVec.2 = ';%-199.752430

% Create posvec from image info
posvec = CreatePosVec( info );

disp( sprintf( '%03d %02d %02d %02d %02d %02d %02d %10.5f %10.5f %10.5f ', ...
    info.InstanceNumber,...
    info.ImageOrientationPatient(4),...
    info.ImageOrientationPatient(5),...
    info.ImageOrientationPatient(6),...
    info.ImageOrientationPatient(1),...
    info.ImageOrientationPatient(2),...
    info.ImageOrientationPatient(3),...
    posvec(1),...
    posvec(2),...
    posvec(3) ) );

header1 = sprintf( '%s\n%s%d\n%s%d\n%s%d\n%s%d\n%s%d\n',...
    s1,s2,anatomicalSliceNumber,s3,info.InstanceNumber-1,s4,sliceNumber,...
    s5,info.Columns,...
    s6,info.Rows );
header2 = sprintf( '%s%f\n%s%f\n%s%f\n%s%f\n',...
     s7,info.SliceThickness,s8,info.ImagePositionPatient(3),s9,info.PixelSpacing(1),s10,info.PixelSpacing(2) );
% tfg - the s8,info.ImagePositionPatient(3) is for axial slice volume data
% we need to change it to (1) for sagittal and (2) for coronal slice data
 
% info.ImagePositionPatient(3) only works for axial slice order
% so correct the line above as we have found that not all DICOM sets will
% have info.SliceLocation set per
% the commented code below
%     s7,info.SliceThickness,s8,info.SliceLocation,s9,info.PixelSpacing(1),s10,info.PixelSpacing(2) );
header3 = sprintf( '%s%f\n%s%f\n%s%f\n%s%f\n%s%f\n%s%f\n%s%f\n',...
    s11,info.ImageOrientationPatient(4),...
    s12,info.ImageOrientationPatient(5),...
    s13,info.ImageOrientationPatient(6),...
    s14,info.ImageOrientationPatient(1),...
    s15,info.ImageOrientationPatient(2),...
    s16,info.ImageOrientationPatient(3),...
    s17,posvec(1),...
    s18,posvec(2),...
    s19,posvec(3) );
header = [ header1 header2 header3 ];
end
