%% Lateral drift correction of Localisation Data based on Fiducial-drift 
%                               TWO clolours

%   This script was develeoped by Timo Rey, LEB, IPHYS-EPFL, Switzerland
%                         Final version December 2019

%   Code is based on Drift_correction-script by Christian Sieben as published in:
%   C. Sieben et al., Nature Methods 15, 777-780 (2018)
%   https://github.com/christian-7
%   This script thus falls under the same license as the original script by C. Sieben [see github].


%% 2-Colour Drift Correction:

clear, clc, close all
% FOVs to do: 1,2,3
FOV = 12;

directory = 'C:\Users\Public\Documents\03_Architecture\02_STORM_2C\03_processedData\';
file_ch1 = '20191212_TR_COS7_BrU-A647_';
file_ch2 = '20191212_TR_COS7_D2-d755_';

% proposed naming structure:
% directory = 'directory of my data\Localizations\';
% file_ch1 = 'date_initials_cellType_target-fluorophore1_';
% file_ch2 = 'date_initials_cellType_target-fluorophore2_';

nameC1 = [directory file_ch1 num2str(FOV) '_1\' file_ch1 num2str(FOV) '_1_MMStack_1_Localizations.csv'];
nameC2 = [directory file_ch2 num2str(FOV) '_1\' file_ch2 num2str(FOV) '_1_MMStack_1_Localizations_affineApplied.csv'];

[filepath_Ch1,name_Ch1,ext_Ch1] = fileparts(nameC1);
[filepath_Ch2,name_Ch2,ext_Ch2] = fileparts(nameC2);

% Load file first channel
cd(filepath_Ch1);locs_Ch1 = dlmread([name_Ch1 ext_Ch1],',',1,0);
cd(filepath_Ch2);locs_Ch2 = dlmread([name_Ch2 ext_Ch2],',',1,0);

locs_Ch1(:,end+1) = 1; % Channel ID, col 9
locs_Ch2(:,end+1) = 2; % Channel ID, col 9

RegionID = size(locs_Ch1,2)+1; % only for Fid variable, col 10

% Load the header and find the right Columns
file = fopen(nameC1);
line = fgetl(file);
h = regexp( line, ',', 'split' );

xCol      = strmatch('x [nm]',h);
yCol      = strmatch('y [nm]',h);
frameCol  = strmatch('frame',h);
deltaXCol = size(locs_Ch1,2)+1; % col 10
deltaYCol = size(locs_Ch1,2)+2; % col 11

display('-- Data loaded --')

%% Combine Channels and select fiducials from Image

allLocs = vertcat(locs_Ch1, locs_Ch2);

pxlsize = 500;

heigth  = round((max(allLocs(:,yCol))-min(allLocs(:,yCol)))/pxlsize);
width   = round((max(allLocs(:,xCol))-min(allLocs(:,xCol)))/pxlsize);
im      = hist3([allLocs(:,xCol),allLocs(:,yCol)],[width heigth]); % heigth x width

% Select rectangles
rect = []; rect2 = [];

figure('Position',[100 50 1000 1000])
f = imagesc(imrotate(im,90),[(max(locs_Ch1(:,frameCol))+max(locs_Ch2(:,frameCol)))*0.6 max(locs_Ch1(:,frameCol))+max(locs_Ch2(:,frameCol))]);
colormap('parula'); colorbar;

while isvalid(f)
  try  rect = getrect;
       rect2 = vertcat(rect2,rect);      
  catch continue
  end

end


fprintf('\n -- Fiducials selected --\n')

% Plot fiducials and average curve rectangles

% Select ROI for both channels
Fid_Ch1 = []; Fid_Ch2 = [];

for i = 1:size(rect2,1);   
xmin = min(allLocs(:,xCol))+ rect2(i,1)*pxlsize;
ymin = max(allLocs(:,yCol)) - rect2(i,2)*pxlsize - (rect2(i,4)*pxlsize) ;
xmax = xmin + (rect2(i,3)* pxlsize);
ymax = ymin + rect2(i,4) * pxlsize;

vx      = find(allLocs(:,xCol)>xmin & allLocs(:,xCol)<xmax);
subset1 = allLocs(vx,1:end);

vy      = find(subset1(:,yCol)>ymin & subset1(:,yCol)<ymax);
subset2 = subset1(vy,1:end);

subset2(:,end+1)=i; % Region ID

Fid_Ch1 = vertcat(Fid_Ch1,subset2(subset2(:,end-1)==1,1:end));
Fid_Ch2 = vertcat(Fid_Ch2,subset2(subset2(:,end-1)==2,1:end));
end

% Plot the fiducials
close all

for i = 1:max(Fid_Ch1(:,end));
    figure('Position',[100 200 400 400])  
    scatter(Fid_Ch1(Fid_Ch1(:,end)==i,frameCol),Fid_Ch1(Fid_Ch1(:,end)==i,yCol),1,'blue');hold on;
    scatter(Fid_Ch2(Fid_Ch2(:,end)==i,frameCol),Fid_Ch2(Fid_Ch2(:,end)==i,yCol),1,'green');
    legend('Ch1','Ch2');
end

%% Select the fiducials and Normalize them to their center of mass
%close all;

selectedFid = [1,6,7,9,10,13,14];

% normalise selected fiducials:
Avg_Ch1x = []; Avg_Ch1y = []; Avg_Ch1 = []; Avg_Ch1frame = [];Avg_Ch1ID = [];

for i = selectedFid;
    target  = find(Fid_Ch1(:,RegionID) == i & Fid_Ch1(:,frameCol)<10000);
    offsetX = median(Fid_Ch1(target,xCol)); offsetY = median(Fid_Ch1(target,yCol)); % median of the first 5000 frames 
    
    Avg_Ch1x        = vertcat(Avg_Ch1x,Fid_Ch1(Fid_Ch1(:,RegionID)==i,xCol)-offsetX);
    Avg_Ch1y        = vertcat(Avg_Ch1y,Fid_Ch1(Fid_Ch1(:,RegionID)==i,yCol)-offsetY);
    Avg_Ch1frame    = vertcat(Avg_Ch1frame,Fid_Ch1(Fid_Ch1(:,RegionID)==i,frameCol));
    Avg_Ch1ID       = vertcat(Avg_Ch1ID,Fid_Ch1(Fid_Ch1(:,RegionID)==i,RegionID)); % Region ID
    
end

Avg_Ch1      = Avg_Ch1x;
Avg_Ch1(:,2) = Avg_Ch1y;
Avg_Ch1(:,3) = Avg_Ch1frame;
Avg_Ch1(:,4) = Avg_Ch1ID; % Region ID

clear Avg_Ch1x Avg_Ch1y Avg_Ch1frame Avg_Ch1ID

Avg_Ch2x = []; Avg_Ch2y = []; Avg_Ch2 = []; Avg_Ch2frame = [];Avg_Ch2ID = [];

for i = selectedFid;
    target = find(Fid_Ch2(:,RegionID)==i & Fid_Ch2(:,frameCol)<10000);
    offsetX = median(Fid_Ch2(target,xCol)); offsetY = median(Fid_Ch2(target,yCol));
    
    Avg_Ch2x        = vertcat(Avg_Ch2x,Fid_Ch2(Fid_Ch2(:,RegionID)==i,xCol)-offsetX);
    Avg_Ch2y        = vertcat(Avg_Ch2y,Fid_Ch2(Fid_Ch2(:,RegionID)==i,yCol)-offsetY);
    Avg_Ch2frame    = vertcat(Avg_Ch2frame,Fid_Ch2(Fid_Ch2(:,RegionID)==i,frameCol));
    Avg_Ch2ID       = vertcat(Avg_Ch2ID,Fid_Ch2(Fid_Ch2(:,RegionID)==i,RegionID)); % Region ID
end

Avg_Ch2      = Avg_Ch2x;
Avg_Ch2(:,2) = Avg_Ch2y;
Avg_Ch2(:,3) = Avg_Ch2frame;
Avg_Ch2(:,4) = Avg_Ch2ID; % Region ID

clear Avg_Ch2x Avg_Ch2y Avg_Ch2frame Avg_Ch2ID

figure('Position', [200 200 400 500])
subplot(2,1,1)
scatter(Avg_Ch1(:,3),Avg_Ch1(:,1),1,'c'); hold on;
scatter(Avg_Ch1(:,3),Avg_Ch1(:,2),1,'r'); hold on;
title('Channel 1'); box on;
legend('X drift', 'Y Drift');

subplot(2,1,2)
scatter(Avg_Ch2(:,3),Avg_Ch2(:,1),1,'c'); hold on;
scatter(Avg_Ch2(:,3),Avg_Ch2(:,2),1,'r'); hold on;
title('Channel 2'); box on;
legend('X drift', 'Y Drift');

%% save metadata:
metaFile = [name_Ch1 '_DC-MetadataFiducials.txt'];
count = 1;

for i = selectedFid;
    fiducials(count,:) = rect2(i,:);
    count = count + 1;
end

dlmwrite(metaFile, fiducials, 'delimiter', ',', 'newline', 'pc');


%% Average the fiducial tracks 

%close all;
frameCol = 3;
xCol     = 1;
yCol     = 2;

% Channel 1

Avg_Ch1_new = []; count = 1;

for i = min(Avg_Ch1(:,frameCol)):max(Avg_Ch1(:,frameCol));      % For all frames
   target = find(Avg_Ch1(:,frameCol) == i);                     % find all fiducials in frame i
   
   if isempty(target);
   else    
   
   Avg_Ch1_new(i,1) = i; % frame
   Avg_Ch1_new(i,2) = mean(Avg_Ch1(target,xCol));               % mean x of all fiducials in frame i
   Avg_Ch1_new(i,3) = mean(Avg_Ch1(target,yCol));               % mean x of all fiducials in frame i
   
   cont = count +1;
   end
end

Avg_Ch1_new(1:min(Avg_Ch1(:,frameCol))-1,:) = [];

% Channel 2

Avg_Ch2_new = []; count = 1;

for i = min(Avg_Ch2(:,frameCol)):max(Avg_Ch2(:,frameCol));      % For all frames
   target = find(Avg_Ch2(:,frameCol) == i);                     % find all fiducials in frame i
   
   if isempty(target);
   else    
   
   Avg_Ch2_new(i,1) = i; % frame
   Avg_Ch2_new(i,2) = mean(Avg_Ch2(target,xCol));               % mean x of all fiducials in frame i
   Avg_Ch2_new(i,3) = mean(Avg_Ch2(target,yCol));               % mean x of all fiducials in frame i
   
   cont = count +1;
   end
end

Avg_Ch2_new(1:min(Avg_Ch2(:,frameCol))-1,:) = [];

figure('Position', [200 200 400 500])
subplot(2,1,1)
scatter(Avg_Ch1_new(:,1),Avg_Ch1_new(:,2),1,'c'); hold on;
scatter(Avg_Ch1_new(:,1),Avg_Ch1_new(:,3),1,'m'); hold on;
title('Channel 1'); box on;
legend('X drift', 'Y Drift');

subplot(2,1,2)
scatter(Avg_Ch2_new(:,1),Avg_Ch2_new(:,2),1,'c'); hold on;
scatter(Avg_Ch2_new(:,1),Avg_Ch2_new(:,3),1,'m'); hold on;
title('Channel 2'); box on;
legend('X drift', 'Y Drift');

%% Spline fit average 
%close all

%%%%%
NbrBins             = 30;    % default = 30
radius              = 90;   % Radius around the fiducial center, default = 250
smoothingFactor     = 100;   % default = 100
startFrame          = [1,100];  % default = 1 <- NOT ALLOWED to be smaller than the first localisation frame.

for allowingLegibility = 1
%%%%%
% Channel 1
    [splineResX,AvgCurveX,pX] = splineFit(Avg_Ch1_new(:,1),Avg_Ch1_new(:,2),NbrBins,radius,smoothingFactor); % (xData,yData,NbrBins,radius,smoothingFactor);
    [splineResY,AvgCurveY,pY] = splineFit(Avg_Ch1_new(:,1),Avg_Ch1_new(:,3),NbrBins,radius,smoothingFactor);

    figure('Position', [200 200 700 500],'NumberTitle', 'off', 'Name', 'Drift correction Ch1')
    subplot(2,3,1)
    scatter(Avg_Ch1_new(:,1),Avg_Ch1_new(:,2),2,'b'), hold on;
    plot(splineResX(:,1),splineResX(:,2),'r.'), hold on;
    axis([0 max(Avg_Ch1_new(:,1)) -radius radius])
    axis square; box on

    subplot(2,3,2)
    plot(AvgCurveX(:,1),AvgCurveX(:,2),'o'); hold on;
    fnplt(csaps(AvgCurveX(:,1),AvgCurveX(:,2),pX/smoothingFactor),'r--')
    legend('noisy data','smoothing spline'), hold off
    axis([0 max(Avg_Ch1_new(:,1)) -radius radius])
    axis square; box on

    subplot(2,3,4)
    scatter(Avg_Ch1_new(:,1),Avg_Ch1_new(:,3),2,'b'), hold on;
    plot(splineResY(:,1),splineResY(:,2),'r.'), hold on;
    axis([0 max(Avg_Ch1_new(:,1)) -radius radius])
    axis square; box on

    subplot(2,3,5)
    plot(AvgCurveY(:,1),AvgCurveY(:,2),'o'); hold on;
    fnplt(csaps(AvgCurveY(:,1),AvgCurveY(:,2),pY/smoothingFactor),'r--')
    legend('noisy data','smoothing spline'), hold off
    axis([0 max(Avg_Ch1_new(:,1)) -radius radius])
    axis square; box on

%%%%%%%%%%%%%%% Correct Channel 1 Averages Tracks

% 1. Calculate the dirft vs. frame curve
    Avg_Ch1_new(:,4) = csaps(AvgCurveX(:,1),AvgCurveX(:,2),pX/100, Avg_Ch1_new(:,1)); % spline fit of the X Ch1 acc frame
    Avg_Ch1_new(:,5) = csaps(AvgCurveY(:,1),AvgCurveY(:,2),pY/100, Avg_Ch1_new(:,1)); % spline fit of the Y Ch1 acc frame

% 2. Subtract the value at startFrame
    Avg_Ch1_new(:,4) = Avg_Ch1_new(:,4)-csaps(AvgCurveX(:,1),AvgCurveX(:,2),pX/100, startFrame(1)); % deltaX
    Avg_Ch1_new(:,5) = Avg_Ch1_new(:,5)-csaps(AvgCurveX(:,1),AvgCurveX(:,2),pX/100, startFrame(1)); % deltaY

% 3. Set everything <startFrame = 0
    Avg_Ch1_new(Avg_Ch1_new(:,1)<startFrame(1),4) = 0; % preallocate delta column to 0
    Avg_Ch1_new(Avg_Ch1_new(:,1)<startFrame(1),5) = 0; % preallocate delta column to 0

%%%%%%%%%%%%%%% Correct Channel 1 Fiducial Tracks

% 1. Calculate the dirft vs. frame curve
    Fid_Ch1(:,deltaXCol+1) = csaps(AvgCurveX(:,1),AvgCurveX(:,2),pX/100, Fid_Ch1(:,frameCol));
    Fid_Ch1(:,deltaYCol+1) = csaps(AvgCurveY(:,1),AvgCurveY(:,2),pY/100, Fid_Ch1(:,frameCol));

% 2. Subtract the value at startFrame
    Fid_Ch1(:,deltaXCol+1) = Fid_Ch1(:,deltaXCol+1)-csaps(AvgCurveX(:,1),AvgCurveX(:,2),pX/100, startFrame(1)); % deltaX
    Fid_Ch1(:,deltaYCol+1) = Fid_Ch1(:,deltaYCol+1)-csaps(AvgCurveY(:,1),AvgCurveY(:,2),pY/100, startFrame(1)); % deltaY

% 3. Set everything <startFrame = 0
    Fid_Ch1(Fid_Ch1(:,frameCol)<startFrame(1),deltaXCol+1) = 0; % preallocate delta column to 0
    Fid_Ch1(Fid_Ch1(:,frameCol)<startFrame(1),deltaYCol+1) = 0;

% 4. Correct the XY coordinates
    Fid_Ch1_DC = [];
    Fid_Ch1_DC(:,xCol)      = Fid_Ch1(:,xCol) - Fid_Ch1(:,deltaXCol+1); % deltaX
    Fid_Ch1_DC(:,yCol)      = Fid_Ch1(:,yCol) - Fid_Ch1(:,deltaYCol+1); % deltaY
    Fid_Ch1_DC(:,frameCol)  = Fid_Ch1(:,frameCol); % deltaY
    Fid_Ch1_DC(:,4)         = Fid_Ch1(:,RegionID); % deltaY

%%%%%%%%%%%%%%% Correct Channel 1 locs

    locs_Ch1_DC = locs_Ch1;

    locs_Ch1_DC(:,deltaXCol) = csaps(AvgCurveX(:,1),AvgCurveX(:,2),pX/100, locs_Ch1_DC(:,frameCol));  % deltaX
    locs_Ch1_DC(:,deltaYCol) = csaps(AvgCurveY(:,1),AvgCurveY(:,2),pY/100, locs_Ch1_DC(:,frameCol));  % deltaY

    locs_Ch1_DC(:,deltaXCol) = locs_Ch1_DC(:,deltaXCol)-csaps(AvgCurveX(:,1),AvgCurveX(:,2),pX/100, startFrame(1)); % deltaX - deltaX(startFrame)
    locs_Ch1_DC(:,deltaYCol) = locs_Ch1_DC(:,deltaYCol)-csaps(AvgCurveY(:,1),AvgCurveY(:,2),pY/100, startFrame(1)); % deltaY - deltaY(startFrame)

    locs_Ch1_DC(locs_Ch1_DC(:,frameCol)<startFrame(1),deltaXCol) = 0; 
    locs_Ch1_DC(locs_Ch1_DC(:,frameCol)<startFrame(1),deltaYCol) = 0;

    locs_Ch1_DC(:,xCol) = locs_Ch1_DC(:,xCol)-locs_Ch1_DC(:,deltaXCol); % substract deltaX from X Col
    locs_Ch1_DC(:,yCol) = locs_Ch1_DC(:,yCol)-locs_Ch1_DC(:,deltaYCol); % substract deltaY from Y Col

    subplot(2,3,3)
    scatter(Avg_Ch1_new(:,1),Avg_Ch1_new(:,2)-Avg_Ch1_new(:,4),1,'b'), hold on;
    axis([0 max(Avg_Ch1_new(:,1)) -radius radius])
    axis square; box on
    title('X trajectory after correction');

    subplot(2,3,6)
    scatter(Avg_Ch1_new(:,1),Avg_Ch1_new(:,3)-Avg_Ch1_new(:,5),1,'b'), hold on;
    axis([0 max(Avg_Ch1_new(:,1)) -radius radius])
    axis square; box on
    title('Y trajectory after correction'); 

    % Channel 2
    [splineResX,AvgCurveX,pX] = splineFit(Avg_Ch2_new(:,1),Avg_Ch2_new(:,2),NbrBins,radius,smoothingFactor);
    [splineResY,AvgCurveY,pY] = splineFit(Avg_Ch2_new(:,1),Avg_Ch2_new(:,3),NbrBins,radius,smoothingFactor);

    figure('Position', [1000 200 700 500],'NumberTitle', 'off', 'Name', 'Drift correction Ch2')
    subplot(2,3,1)
    scatter(Avg_Ch2_new(:,1),Avg_Ch2_new(:,2),2,'b'), hold on;
    plot(splineResX(:,1),splineResX(:,2),'r.'), hold on;
    axis([0 max(Avg_Ch2_new(:,1)) -radius radius])
    axis square; box on
    title('X trajectory before correction');

    subplot(2,3,2)
    plot(AvgCurveX(:,1),AvgCurveX(:,2),'o'); hold on;
    fnplt(csaps(AvgCurveX(:,1),AvgCurveX(:,2),pX/smoothingFactor),'--')
    legend('noisy data','smoothing spline'), hold off
    axis([0 max(Avg_Ch2_new(:,1)) -radius radius])
    axis square; box on

    subplot(2,3,4)
    scatter(Avg_Ch2_new(:,1),Avg_Ch2_new(:,3),2,'b'), hold on;
    plot(splineResY(:,1),splineResY(:,2),'r.'), hold on;
    axis([0 max(Avg_Ch2_new(:,1)) -radius radius])
    axis square; box on
    title('Y trajectory before correction');

    subplot(2,3,5)
    plot(AvgCurveY(:,1),AvgCurveY(:,2),'o'); hold on;
    fnplt(csaps(AvgCurveY(:,1),AvgCurveY(:,2),pY/smoothingFactor),'--')
    legend('noisy data','smoothing spline'), hold off
    axis([0 max(Avg_Ch2_new(:,1)) -radius radius])
    axis square; box on

%%%%%%%%%%%%%%% Correct Channel 2 Averages Tracks

    % 1. Calculate the dirft vs. frame curve
    Avg_Ch2_new(:,4) = csaps(AvgCurveX(:,1),AvgCurveX(:,2),pX/100, Avg_Ch2_new(:,1)); % spline fit of the X Ch1 acc frame
    Avg_Ch2_new(:,5) = csaps(AvgCurveY(:,1),AvgCurveY(:,2),pY/100, Avg_Ch2_new(:,1)); % spline fit of the Y Ch1 acc frame

    % 2. Subtract the value at startFrame
    Avg_Ch2_new(:,4) = Avg_Ch2_new(:,4)-csaps(AvgCurveX(:,1),AvgCurveX(:,2),pX/100, startFrame(2)); % deltaX
    Avg_Ch2_new(:,5) = Avg_Ch2_new(:,5)-csaps(AvgCurveX(:,1),AvgCurveX(:,2),pX/100, startFrame(2)); % deltaY

    % 3. Set everything <startFrame = 0
    Avg_Ch2_new(Avg_Ch2_new(:,1)<startFrame(2),4) = 0; % preallocate delta column to 0
    Avg_Ch2_new(Avg_Ch2_new(:,1)<startFrame(2),5) = 0; % preallocate delta column to 0

%%%%%%%%%%%%%%% Correct Channel 1 Fiducial Tracks

    % 1. Calculate the dirft vs. frame curve
    Fid_Ch2(:,deltaXCol+1) = csaps(AvgCurveX(:,1),AvgCurveX(:,2),pX/100, Fid_Ch2(:,frameCol));
    Fid_Ch2(:,deltaYCol+1) = csaps(AvgCurveY(:,1),AvgCurveY(:,2),pY/100, Fid_Ch2(:,frameCol));

    % 2. Subtract the value at startFrame
    Fid_Ch2(:,deltaXCol+1) = Fid_Ch2(:,deltaXCol+1)-csaps(AvgCurveX(:,1),AvgCurveX(:,2),pX/100, startFrame(2)); % deltaX
    Fid_Ch2(:,deltaYCol+1) = Fid_Ch2(:,deltaYCol+1)-csaps(AvgCurveY(:,1),AvgCurveY(:,2),pY/100, startFrame(2)); % deltaY

    % 3. Set everything <startFrame = 0
    Fid_Ch2(Fid_Ch2(:,frameCol)<startFrame(2),deltaXCol+1) = 0; % preallocate delta column to 0
    Fid_Ch2(Fid_Ch2(:,frameCol)<startFrame(2),deltaYCol+1) = 0;

    % 4. Correct the XY coordinates
    Fid_Ch2_DC = [];
    Fid_Ch2_DC(:,xCol)      = Fid_Ch2(:,xCol) - Fid_Ch2(:,deltaXCol+1); % deltaX
    Fid_Ch2_DC(:,yCol)      = Fid_Ch2(:,yCol) - Fid_Ch2(:,deltaYCol+1); % deltaY
    Fid_Ch2_DC(:,frameCol)  = Fid_Ch2(:,frameCol); % deltaY
    Fid_Ch2_DC(:,4)         = Fid_Ch2(:,RegionID); % deltaY

    %%%%%%%%%%%%%%% Correct Channel 2 locs

    locs_Ch2_DC = locs_Ch2;

    locs_Ch2_DC(:,deltaXCol) = csaps(AvgCurveX(:,1),AvgCurveX(:,2),pX/100, locs_Ch2_DC(:,frameCol));  % deltaX
    locs_Ch2_DC(:,deltaYCol) = csaps(AvgCurveY(:,1),AvgCurveY(:,2),pY/100, locs_Ch2_DC(:,frameCol));  % deltaY

    locs_Ch2_DC(:,deltaXCol) = locs_Ch2_DC(:,deltaXCol)-csaps(AvgCurveX(:,1),AvgCurveX(:,2),pX/100, startFrame(2)); % deltaX - deltaX(startFrame)
    locs_Ch2_DC(:,deltaYCol) = locs_Ch2_DC(:,deltaYCol)-csaps(AvgCurveY(:,1),AvgCurveY(:,2),pY/100, startFrame(2)); % deltaY - deltaY(startFrame)

    locs_Ch2_DC(locs_Ch2_DC(:,frameCol)<startFrame(2),deltaXCol) = 0; 
    locs_Ch2_DC(locs_Ch2_DC(:,frameCol)<startFrame(2),deltaYCol) = 0;

    locs_Ch2_DC(:,xCol) = locs_Ch2_DC(:,xCol)-locs_Ch2_DC(:,deltaXCol); % substract deltaX from X Col
    locs_Ch2_DC(:,yCol) = locs_Ch2_DC(:,yCol)-locs_Ch2_DC(:,deltaYCol); % substract deltaY from Y Col

    subplot(2,3,3)
    scatter(Avg_Ch2_new(:,1),Avg_Ch2_new(:,2)-Avg_Ch2_new(:,4),1,'b'), hold on;
    axis([0 max(Avg_Ch2_new(:,1)) -radius radius])
    axis square; box on
    title('X trajectory after correction');

    subplot(2,3,6)
    scatter(Avg_Ch2_new(:,1),Avg_Ch2_new(:,3)-Avg_Ch2_new(:,5),1,'b'), hold on;
    axis([0 max(Avg_Ch2_new(:,1)) -radius radius])
    axis square; box on
    title('Y trajectory after correction');
end

%% save metafile 2:
metaFile2 = [name_Ch1 '_DC-MetadataFiducials.mat'];
save(metaFile2, 'selectedFid');
save(metaFile2, 'Fid_Ch1', '-append');
save(metaFile2, 'Fid_Ch2', '-append');
save(metaFile2, 'startFrame', '-append');
save(metaFile2, 'smoothingFactor', '-append');
save(metaFile2, 'radius', '-append');
save(metaFile2, 'NbrBins', '-append');

disp('metadata 2 saved');

%% Find CoM of Fiducials (again) and substract from Ch2
%close all;
center_Ch1 = [];center_Ch2 = [];

RegionID = 4;
for i = selectedFid;
    center_Ch1(i+1,1) = median(Fid_Ch1_DC(Fid_Ch1_DC(:,RegionID)==i,xCol));
    center_Ch1(i+1,2) = median(Fid_Ch1_DC(Fid_Ch1_DC(:,RegionID)==i,yCol));
end

for i = selectedFid;
    center_Ch2(i+1,1) = median(Fid_Ch2_DC(Fid_Ch2_DC(:,RegionID)==i,xCol));
    center_Ch2(i+1,2) = median(Fid_Ch2_DC(Fid_Ch2_DC(:,RegionID)==i,yCol));
end

% Delete NaNs
center_Ch1_noNan = [];
center_Ch2_noNan = [];

for i = 1:size(center_Ch1,1);
    if isnan(center_Ch1(i,1))==1;
        
    else
    center_Ch1_noNan(i,1) = center_Ch1(i,1);
    center_Ch1_noNan(i,2) = center_Ch1(i,2);
    end
end

center_Ch1 = [];
center_Ch1(:,1) = nonzeros(center_Ch1_noNan(:,1));
center_Ch1(:,2) = nonzeros(center_Ch1_noNan(:,2));

for i = 1:size(center_Ch2,1);
    if isnan(center_Ch2(i,1))==1;
        
    else
    center_Ch2_noNan(i,1) = center_Ch2(i,1);
    center_Ch2_noNan(i,2) = center_Ch2(i,2);
    end
end

center_Ch2 = [];
center_Ch2(:,1) = nonzeros(center_Ch2_noNan(:,1));
center_Ch2(:,2) = nonzeros(center_Ch2_noNan(:,2));

figure('Position',[500 500 300 300])
scatter(Fid_Ch1(:,xCol),Fid_Ch1(:,yCol),10,'g','filled');hold on;
scatter(Fid_Ch2(:,xCol),Fid_Ch2(:,yCol),10,'r','filled');
scatter(center_Ch1(:,1),center_Ch1(:,2),20,'bo');hold on;
scatter(center_Ch2(:,1),center_Ch2(:,2),20,'bx');hold on;
box on; axis equal;
title('Indentified Fiducial Centers');

figure('Position',[600 300 500 500])
scatter(Fid_Ch1_DC(:,xCol),Fid_Ch1_DC(:,yCol),10,'g','filled');hold on;
scatter(Fid_Ch2_DC(:,xCol),Fid_Ch2_DC(:,yCol),10,'m','filled');
scatter(center_Ch1(:,1),center_Ch1(:,2),20,'bo');hold on;
scatter(center_Ch2(:,1),center_Ch2(:,2),20,'bx');hold on;
box on; axis equal;
title('Indentified Fiducials with DC');

fprintf('\n -- CoM identified --\n')

%% Extract linear Transformation

deltaXY = [];
for i = 1:size(center_Ch1,1);
    deltaXY(i,1) = center_Ch1(i,1) - center_Ch2(i,1);
    deltaXY(i,2) = center_Ch1(i,2) - center_Ch2(i,2);
end

center_Ch2_corr = [];
center_Ch2_corr(:,1) = center_Ch2(:,1) + mean(deltaXY(:,1));
center_Ch2_corr(:,2) = center_Ch2(:,2) + mean(deltaXY(:,2));

TRE = [];
for i = 1:size(center_Ch1,1); 
TRE(:,i) = sqrt((center_Ch1(i,1)-center_Ch2_corr(i,1))^2 + (center_Ch1(i,2)-center_Ch2_corr(i,2))^2);
end

figure
scatter(center_Ch1(:,1),center_Ch1(:,2),10,'g','filled');hold on;
scatter(center_Ch2(:,1) + mean(deltaXY(:,1)),center_Ch2(:,2) + mean(deltaXY(:,2)),10,'m','filled');  
box on; axis equal;
title(['Fid After 2nd trans TRE = ', num2str(mean(TRE)), ' +/- ', num2str(std(TRE)), ' [nm]']);

fprintf('\n -- Linear Transformation extracted --\n')

%% Apply Linear Translation and Save Data 
NameCorrected = [name_Ch1 '_DC' ext_Ch1];

cd(filepath_Ch1);
fileID = fopen(NameCorrected,'w');
fprintf(fileID,[[line,',Channel ID, dx [nm],dy [nm]'] ' \n']);
dlmwrite(NameCorrected,locs_Ch1_DC,'-append');
fclose('all');

fprintf('\n -- Saved Ch1 --\n');
%%
NameCorrected = [name_Ch2 '_DC_corrected' ext_Ch2];

cd(filepath_Ch2);
fileID = fopen(NameCorrected,'w');
fprintf(fileID,[[line,',Channel ID, dx [nm],dy [nm]'] ' \n']);
dlmwrite(NameCorrected,locs_Ch2_DC,'-append');
fclose('all');

fprintf('\n -- Saved Ch2 --\n');

% save remaining deltas for all fiducials
save('todays_deltaXY.mat', 'deltaXY', 'center_Ch2_corr', 'center_Ch1');

fprintf(['\n -- Finished Processing FOV ' num2str(FOV) ' --\n']);
