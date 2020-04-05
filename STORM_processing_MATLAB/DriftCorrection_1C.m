%% Lateral drift correction of Localisation Data based on Fiducial-drift 
%                               ONE clolour

%   This script was develeoped by Timo Rey, LEB, IPHYS-EPFL, Switzerland
%                         Final version December 2019

%   Code is based on Drift_correction-script by Christian Sieben as published in:
%   C. Sieben et al., Nature Methods 15, 777-780 (2018)
%   https://github.com/christian-7
%   This script thus falls under the same licensing as the original code published by C. Sieben [see github].


%% 1-Colour Drift Correction:

clear, clc, close all

FOV = 1;

Path     = 'A:\User\Folder\';
LocDir   = 'Localizations\';
Name     = 'FileName_';
fileName = [Path LocDir Name num2str(FOV) '_1\' Name num2str(FOV) '_1_MMStack_1_Localizations.csv']

[filepath_Ch1,name_Ch1,ext_Ch1] = fileparts(fileName);

% Load file
cd(filepath_Ch1);locs_Ch1 = dlmread([name_Ch1 ext_Ch1],',',1,0);

locs_Ch1(:,end+1) = 1; % Channel ID, col 9

RegionID = size(locs_Ch1,2)+1; % only for Fid variable, col 10

% Load the header and find the right Columns
file = fopen(fileName);
line = fgetl(file);
h = regexp( line, ',', 'split' );

xCol      = strmatch('x [nm]',h);
yCol      = strmatch('y [nm]',h);
frameCol  = strmatch('frame',h);
deltaXCol = size(locs_Ch1,2)+1; % col 10
deltaYCol = size(locs_Ch1,2)+2; % col 11

display('-- Data loaded --')

%% Select fiducials from Image

allLocs = vertcat(locs_Ch1);

pxlsize = 500;

heigth  = round((max(allLocs(:,yCol))-min(allLocs(:,yCol)))/pxlsize);
width   = round((max(allLocs(:,xCol))-min(allLocs(:,xCol)))/pxlsize);
im      = hist3([allLocs(:,xCol),allLocs(:,yCol)],[width heigth]); % heigth x width

% Select rectangles

rect = []; rect2 = [];

figure('Position',[100 100 900 900])
f = imagesc(imrotate(im,90),[(max(locs_Ch1(:,frameCol)))*0.6 max(locs_Ch1(:,frameCol))]);
colormap('parula'); colorbar;

while isvalid(f)

  try  rect = getrect;

       rect2 = vertcat(rect2,rect); 
       
  catch continue
  end

end


fprintf('\n -- Fiducials selected --\n')

% Plot fiducials and average curve rectangles

Fid_Ch1 = [];

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

end

close all

for i = 1:max(Fid_Ch1(:,end));
    
    figure('Position',[100 200 400 400])  
    scatter(Fid_Ch1(Fid_Ch1(:,end)==i,frameCol),Fid_Ch1(Fid_Ch1(:,end)==i,yCol),1,'green');hold on;
    legend('Ch1');%,'Ch2');
    
end

%% Select the fiducials and Normalize them to their center of mass

close all;

selectedFid = [1,2,3];

% save metadata:
metaFile = [name_Ch1 '_DC-Metadata.txt'];
count = 1;

for i = selectedFid;
    fiducials(count,:) = rect2(i,:);
    count = count + 1;
end

dlmwrite(metaFile, fiducials, 'delimiter', ',', 'newline', 'pc');

% normalise selected fiducials:
Avg_Ch1x = []; Avg_Ch1y = []; Avg_Ch1 = []; Avg_Ch1frame = [];Avg_Ch1ID = [];

for i = selectedFid;
    
    target  = find(Fid_Ch1(:,RegionID)==i & Fid_Ch1(:,frameCol)<50000);
    offsetX = median(Fid_Ch1(target,xCol)); offsetY = median(Fid_Ch1(target,yCol)); % median of the first 1000 frames 
    
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

figure('Position', [200 200 400 500])
subplot(2,1,1)
scatter(Avg_Ch1(:,3),Avg_Ch1(:,1),1,'g'); hold on;
scatter(Avg_Ch1(:,3),Avg_Ch1(:,2),1,'r'); hold on;
title('Channel 1'); box on;
legend('X drift', 'Y Drift');

display('-- fiducials normalised --');

%% Average the fiducial tracks - this creates 1 average track from all tracks selected

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


figure('Position', [200 200 700 700])
subplot(2,1,1)
scatter(Avg_Ch1_new(:,1),Avg_Ch1_new(:,2),1,'g'); hold on;
scatter(Avg_Ch1_new(:,1),Avg_Ch1_new(:,3),1,'r'); hold on;
title('Channel 1'); box on;
legend('X drift', 'Y Drift');

display('-- fiducials averaged --');


%% Spline fit average 

NbrBins             = 30;  % default = 30
radius              = 150; % Radius around the fiducial center, default = 150
smoothingFactor     = 100; % default = 100
startFrame          = 1; % default = 1 <- NOT ALLOWED to be smaller than the first localisation frame.

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

Avg_Ch1_new(:,4) = Avg_Ch1_new(:,4)-csaps(AvgCurveX(:,1),AvgCurveX(:,2),pX/100, startFrame); % deltaX
Avg_Ch1_new(:,5) = Avg_Ch1_new(:,5)-csaps(AvgCurveX(:,1),AvgCurveX(:,2),pX/100, startFrame); % deltaY

% 3. Set everything <startFrame = 0

Avg_Ch1_new(Avg_Ch1_new(:,1)<startFrame,4) = 0; % preallocate delta column to 0
Avg_Ch1_new(Avg_Ch1_new(:,1)<startFrame,5) = 0; % preallocate delta column to 0 -> frames before starting frame will not be _DCed

%%%%%%%%%%%%%%% Correct Channel 1 Fiducial Tracks

% 1. Calculate the dirft vs. frame curve

Fid_Ch1(:,deltaXCol+1) = csaps(AvgCurveX(:,1),AvgCurveX(:,2),pX/100, Fid_Ch1(:,frameCol));
Fid_Ch1(:,deltaYCol+1) = csaps(AvgCurveY(:,1),AvgCurveY(:,2),pY/100, Fid_Ch1(:,frameCol));

% 2. Subtract the value at startFrame

Fid_Ch1(:,deltaXCol+1) = Fid_Ch1(:,deltaXCol+1)-csaps(AvgCurveX(:,1),AvgCurveX(:,2),pX/100, startFrame); % deltaX
Fid_Ch1(:,deltaYCol+1) = Fid_Ch1(:,deltaYCol+1)-csaps(AvgCurveY(:,1),AvgCurveY(:,2),pY/100, startFrame); % deltaY

% 3. Set everything <startFrame = 0

Fid_Ch1(Fid_Ch1(:,frameCol)<startFrame,deltaXCol+1) = 0; % preallocate delta column to 0
Fid_Ch1(Fid_Ch1(:,frameCol)<startFrame,deltaYCol+1) = 0;

% 4. Correct the XY coordinates

Fid_Ch1_DC = [];
Fid_Ch1_DC(:,xCol)      = Fid_Ch1(:,xCol) - Fid_Ch1(:,deltaXCol+1); % deltaX
Fid_Ch1_DC(:,yCol)      = Fid_Ch1(:,yCol) - Fid_Ch1(:,deltaYCol+1); % deltaY
Fid_Ch1_DC(:,frameCol)  = Fid_Ch1(:,frameCol); % deltaY
Fid_Ch1_DC(:,4)         = Fid_Ch1(:,RegionID); % deltaY

% Test it
% figure
% scatter(Fid_Ch1_DC(:,frameCol),Fid_Ch1_DC(:,xCol),1,'k');hold on;
% scatter(Fid_Ch1_DC(:,frameCol),Fid_Ch1_DC(:,yCol),1,'r');

%%%%%%%%%%%%%%% Correct Channel 1 locs

locs_Ch1_DC = locs_Ch1;

locs_Ch1_DC(:,deltaXCol) = csaps(AvgCurveX(:,1),AvgCurveX(:,2),pX/100, locs_Ch1_DC(:,frameCol));  % deltaX
locs_Ch1_DC(:,deltaYCol) = csaps(AvgCurveY(:,1),AvgCurveY(:,2),pY/100, locs_Ch1_DC(:,frameCol));  % deltaY

locs_Ch1_DC(:,deltaXCol) = locs_Ch1_DC(:,deltaXCol)-csaps(AvgCurveX(:,1),AvgCurveX(:,2),pX/100, startFrame); % deltaX - deltaX(startFrame)
locs_Ch1_DC(:,deltaYCol) = locs_Ch1_DC(:,deltaYCol)-csaps(AvgCurveY(:,1),AvgCurveY(:,2),pY/100, startFrame); % deltaY - deltaY(startFrame)

locs_Ch1_DC(locs_Ch1_DC(:,frameCol)<startFrame,deltaXCol) = 0; 
locs_Ch1_DC(locs_Ch1_DC(:,frameCol)<startFrame,deltaYCol) = 0;

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

%% Find CoM of Fiducials and filter out NaNs

center_Ch1 = [];

RegionID = 4;
for i = selectedFid;
    
    center_Ch1(i+1,1) = median(Fid_Ch1_DC(Fid_Ch1_DC(:,RegionID)==i,xCol));
    center_Ch1(i+1,2) = median(Fid_Ch1_DC(Fid_Ch1_DC(:,RegionID)==i,yCol));
    
end

% Delete NaNs
center_Ch1_noNan = [];

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

figure('Position',[200 200 300 300])
scatter(Fid_Ch1(:,xCol),Fid_Ch1(:,yCol),10,'g','filled');hold on;
scatter(center_Ch1(:,1),center_Ch1(:,2),20,'bo');hold on;
box on; axis equal;
title('Indentified Fiducial Centers');

fprintf('\n -- CoM identified --\n')

%% Extract linear Transformation

deltaXY = [];
for i = 1:size(center_Ch1,1);
    
    deltaXY(i,1) = center_Ch1(i,1);
    deltaXY(i,2) = center_Ch1(i,2);
    
end

TRE = [];
for i = 1:size(center_Ch1,1);
    
TRE(:,i) = sqrt((center_Ch1(i,1))^2 + (center_Ch1(i,2))^2);

end

figure
scatter(center_Ch1(:,1),center_Ch1(:,2),10,'g','filled');hold on;
box on; axis equal;
title(['Fid After 2nd trans TRE = ', num2str(mean(TRE))]);

fprintf('\n -- Linear Transformation extracted --\n')

%% Apply Linear Translation and Save Data 
%close all;

NameCorrected = [name_Ch1 '_DC' ext_Ch1];

cd(filepath_Ch1);
fileID = fopen(NameCorrected,'w');
fprintf(fileID,[[line,',Channel ID, dx [nm],dy [nm]'] ' \n']);
dlmwrite(NameCorrected,locs_Ch1_DC,'-append');
fclose('all');

fprintf('\n -- Saved Ch1 --\n');

fprintf(['\n -- Finished Processing FOV ' num2str(FOV) ' --\n']);
