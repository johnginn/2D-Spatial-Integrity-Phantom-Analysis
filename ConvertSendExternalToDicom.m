function [volume] = ConvertSendExternalToDicom( filename, someDicomFileName, displayModulo, maxSlice, outputDirectory, seriesNumber,saveDicom,...
    fileDirectory,currentDirectory)
% To view a time series used in gating.
% filename is the file to read
cd(fileDirectory)
fid = fopen( filename, 'rb' );
numberOfImages = 0;
cd(currentDirectory)

fprintf( '%s\n', filename );

for itime = 1:1:maxSlice
    
    % How big is the header.
    % This is many ascii strings of dicom data.
    % We will parse for rows and columns
    headerSize = fread( fid, 1, 'int32' );
    
    % File ends when they write a zero header size
    if headerSize == 0, break, end;
    % Or the header comes back with nothing in it.
    if isempty(headerSize), return, end;
    
    % How big is the data.
    dataSize = fread( fid, 1, 'int32' );

    % Read the ascii header data
    header = fread( fid, headerSize, 'int8' );
    [asciiDicomTags, count] = sscanf( char(header), '%s' );
    
    % Find number of rows in the header
    rowStringLocation = strfind( asciiDicomTags, 'DICOM.NoOfRows' );
    rowStringLength = 15; % DICOM.NoOfRols=
    firstChar = rowStringLocation + rowStringLength;
    lastChar = firstChar + 5;
    if lastChar > length( asciiDicomTags )
        lastChar = length( asciiDicomTags );
    end
    ysizeChar = asciiDicomTags( firstChar: lastChar );
    ysize = sscanf( ysizeChar, '%d' );
    
    % Find number of columns in the header
    colStringLocation = strfind( asciiDicomTags, 'DICOM.NoOfCols' );
    colStringLength = 15; % DICOM.NoOfCols=
    firstChar = colStringLocation + colStringLength;
    lastChar = firstChar + 5;
    if lastChar > length( asciiDicomTags )
        lastChar = length( asciiDicomTags );
    end    
    xsizeChar = asciiDicomTags( firstChar: lastChar );
    xsize = sscanf( xsizeChar, '%d' );
    
    % DICOM.PosVec.0
    colStringLocation = strfind( asciiDicomTags, 'DICOM.PosVec.0' );
    colStringLength = 15; % length(DICOM.PosVec.0=)
    firstChar = colStringLocation + colStringLength;
    lastChar = firstChar + 5;
    if lastChar > length( asciiDicomTags )
        lastChar = length( asciiDicomTags );
    end    
    xsizeChar = asciiDicomTags( firstChar: lastChar );
    xloc(itime) = sscanf( xsizeChar, '%f' );

    % DICOM.PosVec.1
    colStringLocation = strfind( asciiDicomTags, 'DICOM.PosVec.1' );
    colStringLength = 15; % length(DICOM.PosVec.1=)
    firstChar = colStringLocation + colStringLength;
    lastChar = firstChar + 5;
    if lastChar > length( asciiDicomTags )
        lastChar = length( asciiDicomTags );
    end    
    xsizeChar = asciiDicomTags( firstChar: lastChar );
    yloc(itime) = sscanf( xsizeChar, '%f' );

    % DICOM.PosVec.2
    colStringLocation = strfind( asciiDicomTags, 'DICOM.PosVec.2' );
    colStringLength = 15; % length(DICOM.PosVec.2=)
    firstChar = colStringLocation + colStringLength;
    lastChar = firstChar + 5;
    if lastChar > length( asciiDicomTags )
        lastChar = length( asciiDicomTags );
    end    
    xsizeChar = asciiDicomTags( firstChar: lastChar );
    zloc(itime) = sscanf( xsizeChar, '%f' );
    
    % DICOM.SliceThickness = 5.000000 
    colStringLocation = strfind( asciiDicomTags, 'DICOM.SliceThickness' );
    colStringLength = 21; 
    firstChar = colStringLocation + colStringLength;
    lastChar = firstChar + 5;
    if lastChar > length( asciiDicomTags )
        lastChar = length( asciiDicomTags );
    end    
    xsizeChar = asciiDicomTags( firstChar: lastChar );
    sliceThick(itime) = sscanf( xsizeChar, '%f' );

    % DICOM.SliceLocation = 0.000000
    colStringLocation = strfind( asciiDicomTags, 'DICOM.SliceLocation' );
    colStringLength = 20; 
    firstChar = colStringLocation + colStringLength;
    lastChar = firstChar + 5;
    if lastChar > length( asciiDicomTags )
        lastChar = length( asciiDicomTags );
    end    
    xsizeChar = asciiDicomTags( firstChar: lastChar );
    sliceLocation(itime) = sscanf( xsizeChar, '%f' );

    % DICOM.PixelSpacing.0 = 3.515625
    colStringLocation = strfind( asciiDicomTags, 'DICOM.PixelSpacing.0' );
    colStringLength = 21; 
    firstChar = colStringLocation + colStringLength;
    lastChar = firstChar + 5;
    if lastChar > length( asciiDicomTags )
        lastChar = length( asciiDicomTags );
    end    
    xsizeChar = asciiDicomTags( firstChar: lastChar );
    pixelSpacing0(itime) = sscanf( xsizeChar, '%f' );
    
    % DICOM.PixelSpacing.1 = 3.515625
    colStringLocation = strfind( asciiDicomTags, 'DICOM.PixelSpacing.1' );
    colStringLength = 21; 
    firstChar = colStringLocation + colStringLength;
    lastChar = firstChar + 5;
    if lastChar > length( asciiDicomTags )
        lastChar = length( asciiDicomTags );
    end    
    xsizeChar = asciiDicomTags( firstChar: lastChar );
    pixelSpacing1(itime) = sscanf( xsizeChar, '%f' );
    
    fprintf( 'time=%05d position = %e %e %e %e %e %e\n', ...
        itime, xloc(itime), yloc(itime), zloc(itime),...
        sliceLocation(itime), pixelSpacing0(itime), pixelSpacing1(itime) );
    
    % If this is the first image, we must allocate the volume
    if numberOfImages == 0
        volume( 1:ysize, 1:xsize, 1 ) = 0;
        numberOfImages = numberOfImages + 1;
    else
        numberOfImages = numberOfImages + 1;
    end
    
    image = fread( fid, xsize*ysize, 'int16' );
    image = reshape( image, xsize, ysize );
    image = permute( image, [ 2 1 ] );
    
    volume(:,:,numberOfImages) = image;
    
    % Display input image every so often.
    mod100 = mod( numberOfImages, displayModulo );
    
    if mod100 == 0
                
        imagesc( image );
        title( ['time=' num2str(itime) ] );
        colormap gray;
        axis image;
        pause(0.1);
        
    end
            
end

fclose(fid);

% Now write the volume as dicom images.

% Open a dicom image as a reference
try
    image = dicomread( someDicomFileName );
    info = dicominfo( someDicomFileName );
catch ME
    
    fprintf( 'Exception reading dicom file %s\n', someDicomFileName );
    fprintf( '%s\n', ME.message );
    return;
    
end    

% output series will have a new uid
seriesuid = dicomuid;

% output study will also have a new uid.
studyuid = dicomuid;

% dimensions
nsize = size(volume);

% Reset some of the dicom fields for this patient
info.PatientName.FamilyName = 'RogerNana';
info.PatientName.GivenName = 'ConvertedFromSendExternal';
info.PatientID=sprintf( 'VerificationTesting' );
info.PatientSex='M';
info.PatientBirthDate='19700901';
info.Rows = nsize(2);
info.Columns = nsize(1);
info.PixelSpacing = [ pixelSpacing0(1)*10, pixelSpacing1(1)*10 ]; % units are mm
info.SeriesNumber = seriesNumber;
info.StudyInstanceUID = studyuid;
info.SeriesInstanceUID = seriesuid;

% Edited by JGinn
if saveDicom == 1;
    % Now create each slice, changing dicom header each time.
    for islice = 1:1:nsize(3)
        image = volume(:,:,islice);
        % Set the other dicom fields that have to be set for each slice
        info.SliceLocation = sliceLocation(islice);
        info.SliceLocation = info.SliceLocation * 10; % mm
        info.InstanceNumber = islice;
        info.ImagePositionPatient = [ xloc(islice) yloc(islice) zloc(islice) ];
        simage = int16(image);
        
        % Generate a filename and write it to the disk as dicom.
        try
            
            if islice == 1
                mkdir( outputDirectory );
            end
            
            fname = sprintf( '%s\\%s %d %d.dcm' ,...
                outputDirectory, info.PatientID, seriesNumber, islice );
            dicomwrite( simage, fname, info );
            
            fprintf( '%d of %d Wrote %s\n', islice, nsize(3), fname );
            
        catch ME
            
            fprintf( 'Exception writing dicom file %s\n', fname );
            fprintf( '%s\n', ME.message );
            return;
            
        end
        
    end
end

end