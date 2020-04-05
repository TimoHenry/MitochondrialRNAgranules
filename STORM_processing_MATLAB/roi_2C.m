%% Filtering of Localisations from segmented ROIs for 2-colour STORM data

%   This script was develeoped by Timo Rey, LEB, IPHYS-EPFL, Switzerland
%                         Final version March 2020

%   Code is based on particle_segmentation-script by Christian Sieben as published in:
%   C. Sieben et al., Nature Methods 15, 777-780 (2018)
%   https://github.com/christian-7
%   This script thus falls under the same license as the original script by C. Sieben [see github].

% % Input from roiSegmentation: 
% Drift-corrected localisations
% WF
% ROI-data
% delta-XY
% if previously determined: metadata for WF-adjustments etc.

% % Workflow: 
% Load widefield and ROI-data
% Adjust segmentation parameters for WF and turn into binary               needed for new data-sets.
% Save metadata
% Load localisation data:
% Extract rough particle locations from WF - create ROI
% Correction Factor to stretch WF xy dimensions to Loc dimensions
% Filter neighbours with overlapping boxes
% {Keep only ROIs nearby a Fiducial}                                       NOT used currently
% Build box around each Center and copy locs into variable, Cent
% find center of mass of localisations and re-shift ROIs accordingly
% Add back locs of Fiducials
% Correct localisations according to deltaXY of closest fiducial.
% Visually filter ROIs
% keep or discard:
% Save Locs
% Overlay of kept ROIs only
% Overlay of all ROIs

% % Output: 
% Cent5 - Locs of kept ROIs
% kept ROI-coordinates as .m & .csv for later tracing back of clusters to WF
% Locs overlaid on WF.

%% Read Data
clear, clc, close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FOVs to do: 1,3,5,6,9,10,11,12
i = 12; %for single FOV processing

%%%%%%%%%%%%%%%%% Manual Input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Path            = 'C:\Users\Public\Documents\03_Architecture\02_STORM_2C\01_rawData\';
name_Ch1        = '20200219_TR_COS7_BrU-A647_';
name_Ch2        = '20200219_TR_COS7_D2-d755_';

pxl       = 106;                                                           % Pixel size in nm
filetype  = 2;                                                             % 1 for ThunderStorm, 2 for B-Store

%%%%%%%%%%%%%%%%% Find the data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IM_number = i;

WFpath1          = [Path 'Widefield_images/' name_Ch1 'WF' num2str(IM_number)]; 
WF_name1         = [name_Ch1 'WF' num2str(IM_number) '_MMStack_Pos0.ome.tif'];  

WFpath2         = [Path 'Widefield_images/' name_Ch2 'WF' num2str(IM_number)];
WF_name2        = [name_Ch2 'WF' num2str(IM_number) '_MMStack_Pos0.ome.tif'];          

WF_path = WFpath2; WF_name = WF_name2;
ROI_name    = [name_Ch2 'WF' num2str(IM_number) '_roi_conservative.csv'];    

Locpath1        = [Path 'Localisations/' name_Ch1 num2str(IM_number) '_1'];
locName1        = [name_Ch1 num2str(IM_number) '_1_MMStack_1_Localizations_DC.csv'];

Locpath2        = [Path 'Localisations/' name_Ch2 num2str(IM_number) '_1'];
locName2        = [name_Ch2 num2str(IM_number) '_1_MMStack_1_Localizations_affineApplied_DC_corrected.csv'];

savepath        = [Path 'Analysis/'];
savepath_Images = [savepath 'Output/'];
savename        = [name_Ch1 'FOV_' num2str(IM_number) '_extractedParticles'];

% load delta values:
cd(Locpath2);
load('todays_deltaXY.mat');

fprintf('\n -- Path and File information loaded --\n')

%% Load widefield and ROI-data
cd(WF_path);
ICh = imread(WF_name);                                                     % Widefield

ROIformat = '%*s%*s%*s%f%f%f%f%*s%*s%*s%*s%*s%*s%*s%[^\n\r]';              % ROI-data
RoiID = fopen(ROI_name,'r');
ROIarray = textscan(RoiID, ROIformat, 'Delimiter', ',', 'HeaderLines' ,1, 'ReturnOnError', false);
fclose(RoiID);

fprintf('\n -- WF & ROI data loaded --\n')

% turn ROIs into binary
roiX = ROIarray{:, 2};                                                     % note: pay attention to matlab coordinates
roiY = ROIarray{:, 1};
roiWidth = ROIarray{:, 4};
roiHeight = ROIarray{:, 3};

ROImap = zeros(length(ICh), 'logical'); %create empty table
%
for n = 1:length(roiX);    
    x = 0;
    while x < roiWidth(n);
        y = 0;
        while y < roiHeight(n);
            ROImap(roiX(n) + x, roiY(n) + y) = 1;
            y = y+1;
        end
        x = x + 1;
    end
end
% tidy up (clear variables)
vars = {'x', 'y', 'roiX', 'roiY', 'roiWidth', 'roiHeight', 'n' 'ROIformat', 'RoiID'};
clear(vars{:});
clear vars;

fprintf('\n -- ROIs binarised --\n')

%% Adjust segmentation parameters for WF and turn into binary              needed for new data-sets.
% minWF = 800;
% maxWF = 10000;
% contrast(1,1) = 0.05; contrast(1,2) = 0.25;
% background(1,1) = 0.15; background(1,2) = 0.8;
% binThreshold = 0.25;
load([name_Ch1 'FOV_' num2str(IM_number) '_Segmentation-Metadata.mat']);   % load saved parameters 

% blur the image
I = ICh;
figure('Position',[10 600 500 500],'name','Raw GFP Image'), imshow(I,[minWF maxWF],'InitialMagnification','fit');

% adjust the contrast of the raw image
I2 = imadjust(I,[contrast(1,1) contrast(1,2)],[]);
figure('Position',[600 600 500 500],'name','Image after adjusted contrast'), imshow(I2,'InitialMagnification','fit');

% lowpass filter of size and gaussian blur sigma, [lowpass filter] sigma
G = fspecial('gaussian',[3 3], 100);
imG = imfilter(I2,G,'same');
figure('Position',[1200 600 500 500],'name','Image after gaussian blurring'), imshow(imG,'InitialMagnification','fit');

% adjust the background
I3 = imadjust(imG,[background(1,1) background(1,2)],[]);
figure('Position',[10 10 500 500],'name','Image after adjusted background'), imshow(I3,'InitialMagnification','fit');

% Make binary image
bin = im2bw(I3,binThreshold);
figure('Position',[600 10 500 500],'name','Binary WF-image'),imshow(bin,'InitialMagnification','fit');

% Multiply WF-binary with ROI-binary
bin2 = times(bin,ROImap);
figure('Position',[1200 10 500 500],'name','Final binary'),imshow(bin2,'InitialMagnification','fit');

% show ROIs
figure('Position',[1200 600 500 500],'name','Final binary'),imshow(ROImap,'InitialMagnification','fit');

[B,L,N,A] = bwboundaries(bin2); % B - connectivity

%% save metadata:                                                          in case new parameters were determined
metaFile = [name_Ch1 'FOV_' num2str(IM_number) '_Segmentation-Metadata.mat'];% save(metaFile, 'minWF');
save(metaFile, 'maxWF', '-append');
save(metaFile, 'contrast', '-append');
save(metaFile, 'background', '-append');
save(metaFile, 'binThreshold', '-append');
close all;

%% Load localisation data:
cd(Locpath1);
locs_Ch1=dlmread(locName1,',',1,0);

cd(Locpath2);
locs_Ch2=dlmread(locName2,',',1,0);

% Load header
cd(Locpath1);
file    = fopen(locName1);
line    = fgetl(file);
h       = regexp( line, ',', 'split' );

if filetype == 1;   
    xCol       = strmatch('"x [nm]"',h);
    yCol       = strmatch('"y [nm]"',h);
    LLCol      = strmatch('"loglikelihood"',h);

else
    xCol       = strmatch('x [nm]',h);
    yCol       = strmatch('y [nm]',h);
    LLCol      = strmatch('loglikelihood',h);
end

fprintf('\n -- Localisations loaded --\n')

%% Extract rough particle locations from WF - create ROI
% set ROI-size
box_size = 15;                                                            % Diameter of box in pixels

%Find the center of each ROI and transform into an X,Y coordinate 
Center=[];
for k=1:length(B)
    boundary    = B{k};
	Center(k,1) = (((max(B{k,1}(:,1))-min(B{k,1}(:,1)))/2)+min(B{k,1}(:,1)))*(pxl); % Center of the segmented spot in nm
	Center(k,2) = (((max(B{k,1}(:,2))-min(B{k,1}(:,2)))/2)+min(B{k,1}(:,2)))*(pxl); % Center of the segmented spot in nm
    Center(k,3) = (box_size/2);
end

% Extract the integrated intensity of both WF-channels for each ROI
cd(WFpath1);
ICh1   =   imread(WF_name1);
intI1 = []; Particles_WF1  = {};

cd(WFpath2);
ICh2   =   imread(WF_name2);
intI2 = []; Particles_WF2  = {};

for i = 1:length(B);
    intI1(i,1) = sum(sum(ICh1(min(B{i,1}(:,1)):max(B{i,1}(:,1)),min(B{i,1}(:,2)):max(B{i,1}(:,2)))));
    Particles_WF1{i,1} = ICh1(int16((Center(i,1)./pxl)-Center(i,3)):int16((Center(i,1)./pxl)+Center(i,3)), int16((Center(i,2)./pxl)-Center(i,3)):int16((Center(i,2)./pxl)+Center(i,3)));
    
    intI2(i,1) = sum(sum(ICh2(min(B{i,1}(:,1)):max(B{i,1}(:,1)),min(B{i,1}(:,2)):max(B{i,1}(:,2)))));
    Particles_WF2{i,1} = ICh2(int16((Center(i,1)./pxl)-Center(i,3)):int16((Center(i,1)./pxl)+Center(i,3)), int16((Center(i,2)./pxl)-Center(i,3)):int16((Center(i,2)./pxl)+Center(i,3)));
end
% add box_size to metadata:
save(metaFile, 'box_size', '-append');
disp(' -- boxes built -- ');

%% Correction Factor to stretch WF xy dimensions to Loc dimensions
% Plots of Locs differ from WF bcs of DC & loss of edge-Locs, and are defined by Loc-min & -max.
% CF 'shrinks' the WF to match Loc-plot-dimensions by shifting the ROI centers
center2 = [];
CFX = (max(locs_Ch2(:,xCol)/pxl))./size(I);
CFY = (max(locs_Ch2(:,yCol)/pxl))./size(I);                                

center2(:,1)    =   Center(:,2)*CFX(:,1);                                  % Center of the segmented spot in nm
center2(:,2)    =   Center(:,1)*CFY(:,1);                                  % Center of the segmented spot in nm
center2(:,3)    =   Center(:,3);                                           % width (=diameter) of cluster
Center          =   center2;

fprintf('\n -- %i ROIs found --\n', size(Center,1))
%% Filter neighbours with overlapping boxes
% to avoid clustering particular clusters multiple times.
Distance = []; center2 = [];
count = 1;

for i = 1:size(Center,1);                                                  % for each ROI
    % Check if radii overlap
    otherClusters = Center(setdiff(1:size(Center,1),i),:);                 % find all other spots from the list (Center)
    
    for j = 1:size(otherClusters,1);                                       % for each of these other ROIs
        Distance(j,1) = sqrt((Center(i,1)-otherClusters(j,1))^2+(Center(i,2)-otherClusters(j,2))^2); % find distance between the centers [a^2+b^2=c^2]
    
        if Distance(j,1)>(box_size)*pxl;                                   % if 2 ROI centers have overlapping boxes:
            Distance(j,2) = 0;                                             % assign "boolean-column" 0
        else
            Distance(j,2) = 1;
        end
    end
    
    if sum(Distance(:,2))==0;
        center2(count,1) =  Center(i,1);  
        center2(count,2) =  Center(i,2);  
        center2(count,3) =  Center(i,3);
        
        count = count + 1;
    else end 
end

Center = center2;

fprintf('\n -- Overlapping ROIs filtered: %i ROIs left --\n', size(Center,1))

%% Keep only ROIs nearby a Fiducial
% To correct chromatic shifts, use near-by fiducials. 

% NOT USED CURRENTLY [very large radius to include all ROIs]
% specify fiducial-distance threshold [um]:
radius = 50;
pxl = 106                                                                  % pxl size [nm]
proximity = radius*1000/pxl;                                               % radius in pxls 

% for every ROI-center, find distances to Fiducial-centers (in nm):
FidDistance = pdist2(Center(:,1:2), center_Ch2_corr);

Distance = []; center2 = [];
count = 1;

% Filter ROIs by proximity to Fiducials:
 for i = 1:size(FidDistance,1)                                        % for every ROI
     
    for m = 1:length(FidDistance(i,:))                                     % for each distance to a Fiducial
        if FidDistance(i,m) < proximity*pxl;                               % check if distance is small enough
            Distance(m,1) = 1;                                             % assign "boolean-column" 1
        else 
            Distance(m,1) = 0;                                             % assign "boolean-column" 0
        end
    end
    
    if sum(Distance(:,1)) > 0;                                             % if any distance is small enough, keep the ROI
        center2(count,1) =  Center(i,1);  
        center2(count,2) =  Center(i,2);  
        center2(count,3) =  Center(i,3);

        count = count + 1;
    else end 
 end
 
center3 = Center;
Center = center2;
save(metaFile, 'radius', '-append');

fprintf('\n -- Foci left within Fiducial-radius: %i --\n', length(Center(:,1)))
%% re-initiate to adapt radius-value above
%Center = center3;
%% Build box around each Center and copy locs into variable, Cent
Cent1={}; 
count = 1;

for i = 1:length(Center(:,1));
    disp(i);
    
    vx1 = Center(i,1)+Center(i,3)*pxl;                                     % build box by adding "spot-radius" [note, this was amplified] to center
    vx2 = Center(i,1)-Center(i,3)*pxl;
    vy1 = Center(i,2)+Center(i,3)*pxl;
    vy2 = Center(i,2)-Center(i,3)*pxl;

    target_Ch1 = find(locs_Ch1(:,xCol)>vx2 & locs_Ch1(:,xCol)<vx1 & locs_Ch1(:,yCol)>vy2 & locs_Ch1(:,yCol)<vy1);
    target_Ch2 = find(locs_Ch2(:,xCol)>vx2 & locs_Ch2(:,xCol)<vx1 & locs_Ch2(:,yCol)>vy2 & locs_Ch2(:,yCol)<vy1);

    Cent1{count,1} = locs_Ch1(target_Ch1,1:end);                            % the localisations
    Cent1{count,2} = length(target_Ch1);                                    % number of localisations
    Cent1{count,3} = intI1(i);                                              % WF intensity
    Cent1{count,4} = Particles_WF1{i,1};                                    % locations of the box around the particles in the WF
    Cent1{count,5} = locs_Ch2(target_Ch2,1:end);
    Cent1{count,6} = length(target_Ch2);
    Cent1{count,7} = intI2(i);
    Cent1{count,8} = Particles_WF2{i,1};
      
    count = count + 1;
end 

fprintf('\n -- %i Particles selected from localisation dataset --\n',length(Cent1(:,1)))

%% find center of mass of localisations and re-shift ROIs accordingly
% to adjust ROI-location around clusters of Localisations
CoM = []; Cent={}; count = 1; boxes = [];
segC = 5; % to find locs from channel 2

for i = 1:length(Cent1(:,1));
    disp(i);
    
    CoM(i,1)    =   median(Cent1{i,segC}(:,1));                            % Center of the segmented spot in nm
    CoM(i,2)    =   median(Cent1{i,segC}(:,2));
    CoM(i,3)    =   Center(i,3);                                           % Box 'radius'
    
    % re-build shifted box:
    vx1 = CoM(i,1)+CoM(i,3)*pxl;
    vx2 = CoM(i,1)-CoM(i,3)*pxl;
    vy1 = CoM(i,2)+CoM(i,3)*pxl;
    vy2 = CoM(i,2)-CoM(i,3)*pxl;

    % filter with "spot-boxes":
    target_Ch1 = find(locs_Ch1(:,xCol)>vx2 & locs_Ch1(:,xCol)<vx1 & locs_Ch1(:,yCol)>vy2 & locs_Ch1(:,yCol)<vy1);
    target_Ch2 = find(locs_Ch2(:,xCol)>vx2 & locs_Ch2(:,xCol)<vx1 & locs_Ch2(:,yCol)>vy2 & locs_Ch2(:,yCol)<vy1);
    
    Cent{count,1} = locs_Ch1(target_Ch1,1:end);                            % the localisations
    Cent{count,2} = length(target_Ch1);                                    % number of localisations
    Cent{count,3} = intI1(i);                                              % WF intensity
    Cent{count,4} = Particles_WF1{i,1};                                    % locations of the box around the particles in the WF
    Cent{count,5} = locs_Ch2(target_Ch2,1:end);
    Cent{count,6} = length(target_Ch2);
    Cent{count,7} = intI2(i);
    Cent{count,8} = Particles_WF2{i,1};
    % add ROI-descriptors:
    Cent{count,9} = box_size;
    Cent{count,10} = vx1;
    Cent{count,11} = vy1;
    
    count = count + 1;
end
fprintf('\n -- %i Particles selected after shifting box --\n',length(Cent1(:,1)))

%% Add back locs of Fiducials
Cent4 = Cent;
count = length(Cent(:,1))+1;

for i = 1:length(center_Ch2_corr(:,1));
    disp(i);    
    % build boxes from Fiducial center:
    vx1 = center_Ch2_corr(i,1)+Center(1,3)*pxl;                            % Center(i,3) = Box 'radius'
    vx2 = center_Ch2_corr(i,1)-Center(1,3)*pxl;
    vy1 = center_Ch2_corr(i,2)+Center(1,3)*pxl;
    vy2 = center_Ch2_corr(i,2)-Center(1,3)*pxl;
    
    % filter with "spot-boxes":
    target_Ch1 = find(locs_Ch1(:,xCol)>vx2 & locs_Ch1(:,xCol)<vx1 & locs_Ch1(:,yCol)>vy2 & locs_Ch1(:,yCol)<vy1);
    target_Ch2 = find(locs_Ch2(:,xCol)>vx2 & locs_Ch2(:,xCol)<vx1 & locs_Ch2(:,yCol)>vy2 & locs_Ch2(:,yCol)<vy1);
    
    Cent4{count,1} = locs_Ch1(target_Ch1,1:end);                            % the localisations
    Cent4{count,5} = locs_Ch2(target_Ch2,1:end);

    count = count + 1;    
end

fprintf('\n -- %i Fiducials were added --\n', i)

%% Correct localisations according to deltaXY of closest fiducial.
CoM = []; count = 1; Cent3 = Cent4; 

for i = 1:length(Cent3(:,1));                                              % for every ROI in Cent3
    
    CoM(i,1)    =   median(Cent3{i,5}(:,1));                               % Find center of Mass of Locs in Ch2 (in nm)
    CoM(i,2)    =   median(Cent3{i,5}(:,2));

    FidDistance = pdist2(CoM(i,1:2), center_Ch2_corr);                     % find closest fiducial
    minDist = min(FidDistance);
    myVariable = find(FidDistance == minDist);
    % note: assumes that did not already correct with delta before)
    myX = deltaXY(myVariable(1,1), 1);                                    % calculate remaining Channel-shift in X & Y
    myY = deltaXY(myVariable(1,1), 2);
    
    % correct (shift by myX & myY) Channel 2 localisations
    Cent3{i,5}(:,1) = Cent4{i,5}(:,1)+myX;                                 % the localisations
    Cent3{i,5}(:,2) = Cent4{i,5}(:,2)+myY;
    
end

fprintf('\n -- remaining Channel-shift was corrected. --\n')

%% Visually filter ROIs:
count = 1; 
keptROIs = []; Cent5 = {};
close all

% Plot locs onto WF:
for i = 1:size(Cent,1);                                                    % for every ROI

    % crop WF
    box = 3*box_size;
    xCoord = Cent3{i,10}/pxl-box/2;
    yCoord = Cent3{i,11}/pxl-box/2;
    % display WF
    figure('name', 'Check Loc-shift', 'Position', [500 500 500 500]);
    cropped_img = imcrop(ICh, [xCoord yCoord box box]);
    imshow(cropped_img, [minWF maxWF]); hold on;
    % plot Locs Ch2:
    Cent2=[];
    Cent2(:,1) = Cent3{i,5}(:,xCol);             % divide by correction factor to almost align
    Cent2(:,2) = Cent3{i,5}(:,yCol);             % also correct for estimated shift
    scatter(Cent2(:,1)/pxl-xCoord,Cent2(:,2)/pxl-yCoord,1,'green');
    % Ch1:
    Cent1=[];
    Cent1(:,1) = Cent3{i,1}(:,xCol);             % divide by correction factor to almost align
    Cent1(:,2) = Cent3{i,1}(:,yCol);             % also correct for estimated shift
    scatter(Cent1(:,1)/pxl-xCoord,Cent1(:,2)/pxl-yCoord,1,'blue');
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    %hold on;
end

%%







%% keep or discard:
% specify which ones to keep:
keep = [2,3];


keptROIs = [];
for i = keep
    disp(['ROI ' num2str(i) ' was kept.']);
	% add retained ROIs into new variable:
	keptROIs{count,1} = Cent3{i,10};
	keptROIs{count,2} = Cent3{i,11};
	keptROIs{count,3} = box_size;
	% add retained Locs into new variable:
	Cent5(count,:) = Cent3(i,:);
	count = count + 1;
end

cd(savepath_Images);
% save ROIs:
outFile1 = [name_Ch1 'FOV_' num2str(IM_number) '_ROIsKept.mat'];
save(outFile1, 'keptROIs');
outFile2 = [name_Ch1 'FOV_' num2str(IM_number) '_ROIsKept.csv'];
dlmwrite(outFile2, keptROIs, 'delimiter', ',');

%% Save Locs
cd(savepath_Images);
save(savename,'Cent5','-v7.3');    

fprintf('\n -- %x ROIs were retained for clustering --\n', size(keptROIs, 1));

%% Overlay of kept ROIs only
Cent_wo_duplicates = Cent5;

figure('Position',[150 150 400 400],'name',['Retained ROIs from Ch1 on WF']);
imshow(ICh,[minWF maxWF]); hold on;

for i = 1:size(Cent_wo_duplicates,1);

% for 2nd color:
Cent2=[];
Cent2(:,1) = Cent_wo_duplicates{i,5}(:,xCol)/CFX(:,1);             % divide by correction factor to almost align
Cent2(:,2) = Cent_wo_duplicates{i,5}(:,yCol)/CFY(:,1);             % also correct for estimated shift
scatter(Cent2(:,1)/pxl,Cent2(:,2)/pxl,1,'green');

% same for first color:
Cent1=[];
Cent1(:,1) = Cent_wo_duplicates{i,1}(:,xCol)/CFX(:,1);             % divide by correction factor to almost align
Cent1(:,2) = Cent_wo_duplicates{i,1}(:,yCol)/CFY(:,1);             % also correct for estimated shift
scatter(Cent1(:,1)/pxl,Cent1(:,2)/pxl,1,'blue');

hold on;
end

cd(savepath_Images);
savename_Images = [name_Ch1 'FOV_' num2str(IM_number) '_overlay_kept.fig'];
savefig(savename_Images);

%% Overlay of all ROIs
Cent_wo_duplicates = Cent3;

figure('Position',[150 150 400 400],'name',['Retained ROIs from Ch1 on WF']);
imshow(ICh,[minWF maxWF]); hold on;

for i = 1:size(Cent_wo_duplicates,1);
% for 2nd color:
Cent2=[];
Cent2(:,1) = Cent_wo_duplicates{i,5}(:,xCol)/CFX(:,1);             % divide by correction factor to almost align
Cent2(:,2) = Cent_wo_duplicates{i,5}(:,yCol)/CFY(:,1);             % also correct for estimated shift
scatter(Cent2(:,1)/pxl,Cent2(:,2)/pxl,1,'green');

% same for first color:
Cent1=[];
Cent1(:,1) = Cent_wo_duplicates{i,1}(:,xCol)/CFX(:,1);             % divide by correction factor to almost align
Cent1(:,2) = Cent_wo_duplicates{i,1}(:,yCol)/CFY(:,1);             % also correct for estimated shift
scatter(Cent1(:,1)/pxl,Cent1(:,2)/pxl,1,'blue');

hold on;
end

cd(savepath_Images);
savename_ImagesAll = [name_Ch1 'FOV_' num2str(IM_number) '_overlay_all.fig'];
savefig(savename_ImagesAll);

fprintf('\n -- All Figs saved. cycle back up for more. --\n')
