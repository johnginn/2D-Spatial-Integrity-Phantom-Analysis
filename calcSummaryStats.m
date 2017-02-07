% Function to calculate the summary statistics for each individual slice
% thickness scan.
%
% Input:
% radius (mm) The distance from isocenter
% dev2D (mm) The deviation between the found location and the expected location
% radiusSearch1 (mm) Distance defining the first search region
% radiusSearch2 (mm) Distance defining the second search region
% radiusSearch3 (mm) Distance defining the third search region
% radiusThreshold1 (mm) Deviation tolerance for the first search region
% radiusThreshold2 (mm) Deviation tolerance for the second search region
% radiusThreshold3 (mm) Deviation tolerance for the third search region
% allThreshold (mm) Threshold for all spheres
% saveIndStats y=1,n=0 Whether or not to save the individual stats
% saveFilename The filename for saving the data
% 
% Output:
% avgDevRad1 (mm) The average deviation for cylinders in radiusSearch1
% avgDevRad2 (mm) The average deviation for cylinders in radiusSearch2
% avgDevRad3 (mm) The average deviation for cylinders in radiusSearch3
% maxDevRad1 (mm) The maximum deviation for cylinders in radiusSearch1
% maxDevRad2 (mm) The maximum deviation for cylinders in radiusSearch2
% maxDevRad3 (mm) The maximum deviation for cylinders in radiusSearch3
%
% John Ginn
% Created: 12/5/16
% Modified: 12/5/16

function [avgDevRad1,avgDevRad2,avgDevRad3,maxDevRad1,maxDevRad2,maxDevRad3] = ...
    calcSummaryStats(radius,dev2D,radiusSearch1,radiusSearch2,radiusSearch3,...
    radiusThreshold1,radiusThreshold2,radiusThreshold3,allThreshold,saveIndStats,saveFilename)

countRad1 = 0;
countRad2 = 0;
countRad3 = 0;
countRadPass1 = 0;
countRadPass2 = 0;
countRadPass3 = 0;
countPassAll = 0;
avgDevRad1 = 0;
avgDevRad2 = 0;
avgDevRad3 = 0;
maxDevRad1 = 0;
maxDevRad2 = 0;
maxDevRad3 = 0;
for step = 1:length(radius)
    currRadius = radius(step);
    currDev = dev2D(step);
    if currRadius < radiusSearch1
        countRad1 = countRad1 + 1;
        avgDevRad1 = avgDevRad1 + currDev;
        % check max deviation
        if currDev > maxDevRad1
            maxDevRad1 = currDev;
        end
        if currDev < radiusThreshold1
            countRadPass1 = countRadPass1 + 1;
        end
    end
    if currRadius < radiusSearch2
        countRad2 = countRad2 + 1;
        avgDevRad2 = avgDevRad2 + currDev;
        % check max deviation
        if currDev > maxDevRad2
            maxDevRad2 = currDev;
        end
        if currDev < radiusThreshold2
            countRadPass2 = countRadPass2 + 1;
        end
    end
    if currRadius < radiusSearch3
        countRad3 = countRad3 + 1;
        avgDevRad3 = avgDevRad3 + currDev;
        % check max deviation
        if currDev > maxDevRad3
            maxDevRad3 = currDev;
        end
        if currDev < radiusThreshold3
            countRadPass3 = countRadPass3 + 1;
        end
    end
    if currDev < allThreshold
       countPassAll = countPassAll + 1; 
    end
end

% calculate average deviation
avgDevRad1 = avgDevRad1/countRad1;
avgDevRad2 = avgDevRad2/countRad2;
avgDevRad3 = avgDevRad3/countRad3;
percentPassRad1 = countRadPass1/countRad1*100;
percentPassRad2 = countRadPass2/countRad2*100;
percentPassRad3 = countRadPass3/countRad3*100;
percentPassAll = countPassAll/length(radius)*100;
avgDevAll = sum(dev2D)/length(dev2D);
maxDevAll = max(dev2D);

if saveIndStats == 1
    summaryData = {'Average Deviation for all cylinders',...
        avgDevAll,' (mm)';...
        strcat(['Average Deviation within ',num2str(radiusSearch1),' (mm) of isocenter']),...
        avgDevRad1,' (mm)';...
        strcat(['Average Deviation within ',num2str(radiusSearch2),' (mm) of isocenter']),...
        avgDevRad2,' (mm)';...
        strcat(['Average Deviation within ',num2str(radiusSearch3),' (mm) of isocenter']),...
        avgDevRad3,' (mm)';...
        strcat(['Percent of spheres within ',num2str(radiusSearch1),' (mm) of isocenter w/ less than ',num2str(radiusThreshold1),' (mm) deviation',])...
        percentPassRad1,' (%)';...
        strcat(['Percent of spheres within ',num2str(radiusSearch2),' (mm) of isocenter w/ less than ',num2str(radiusThreshold2),' (mm) deviation',])...
        percentPassRad2,' (%)';...
        strcat(['Percent of spheres within ',num2str(radiusSearch3),' (mm) of isocenter w/ less than ',num2str(radiusThreshold3),' (mm) deviation',])...
        percentPassRad3,' (%)';...
        strcat(['Percent of all spheres w/ less than ',num2str(allThreshold),' (mm) deviation'])...
        percentPassAll,' (%)';...
        'Maximum Deviation for all cylinders',...
        maxDevAll,' (mm)';...
        strcat(['Maximum Deviation within ',num2str(radiusSearch1),' (mm) of isocenter']),...
        maxDevRad1,' (mm)';...
        strcat(['Maximum Deviation within ',num2str(radiusSearch2),' (mm) of isocenter']),...
        maxDevRad2,' (mm)';...
        strcat(['Maximum Deviation within ',num2str(radiusSearch3),' (mm) of isocenter']),...
        maxDevRad3,' (mm)'};
        xlswrite(saveFilename,summaryData);
end

end