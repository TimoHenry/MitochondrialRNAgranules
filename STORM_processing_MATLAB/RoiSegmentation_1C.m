%% Segmentation of ROIs (TS format)

%   This script was develeoped by Timo Rey, LEB, IPHYS-EPFL, Switzerland
%                         Final version December 2019

%   Code is based on ParticleFilter-script by Christian Sieben as published in:
%   C. Sieben et al., Nature Methods 15, 777-780 (2018)
%   https://github.com/christian-7
%   This script thus falls under the same licensing as the original code published by C. Sieben [see github].

% Input: 
%         individual localization files   --> loc path and name
%         individual WF images            --> WF path and name
%         ROIs of in-focus MRGs           --> ROI name (path = WF path)

% Workflow: 
%    Load WF images & ROIs
%    Binarize
%    Load and segment localisations
%    Produce overlay images and save the extraced particles

% Output: 
%       Variable Cent
%    Cent{i,1} = locs_Ch1(target_Ch1,1:9);
%    Cent{i,2} = length(target_Ch1);
%    Cent{i,3} = intI(i);
%    Cent{i,4} = Particles_WF{i,1};


%% Read Data

%for i = [1,2,3];                                                              % for batch processing
clear, clc, close all

%%%%%%%%%%%%%%%%% Manual Input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i          = 1;                                                           % for single FOV processing

Path       = 'A:\User\Folder\';
name_Ch1   = 'FileName';
savepath   = [Path];



%%%%%%%%%%%%%%%%% Find the data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(i);
IM_number = i;
pxl       = 106;                                                           % Pixel size in nm
filetype  = 1;                                                             % 1 for ThunderStorm, 2 for B-Store


WFpath          = [Path 'rawData\' name_Ch1 'WF' num2str(IM_number)]; 
WF_name         = [name_Ch1 'WF' num2str(IM_number) '_MMStack_Pos0.ome.tif'];  

ROI_name    = [name_Ch1 'WF' num2str(IM_number) '_roi.csv'];

Locpath1        = [Path 'Localizations\' name_Ch1 num2str(IM_number) '_1'];
locName1        = [name_Ch1 num2str(IM_number) '_1_MMStack_1_Localizations_DC.csv'];

savepath_Images = [savepath 'Images\'];
savename        = [name_Ch1 'FOV_' num2str(IM_number) '_extractedParticles'];

fprintf('\n -- Path and File information loaded --\n')

%% Load widefield and ROI-data

cd(WFpath);
ICh1       =   imread(WF_name);                                            % Widefield

ROIformat  = '%*s%*s%*s%f%f%f%f%*s%*s%*s%*s%*s%*s%*s%[^\n\r]';             % ROI-data
RoiID      = fopen(ROI_name,'r');
ROIarray   = textscan(RoiID, ROIformat, 'Delimiter', ',', 'HeaderLines' ,1, 'ReturnOnError', false);

fclose(RoiID);

fprintf('\n -- WF & ROI data loaded --\n')

% turn ROIs into binary
roiX       = ROIarray{:, 2};                                               % note: pay attention to matlab coordinates
roiY       = ROIarray{:, 1};
roiWidth   = ROIarray{:, 4};
roiHeight  = ROIarray{:, 3};

ROImap     = zeros(length(ICh1), 'logical');                               %create empty table

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

% tidy up
vars = {'x', 'y', 'roiX', 'roiY', 'roiWidth', 'roiHeight', 'n' 'ROIformat', 'RoiID'};
clear(vars{:});
clear vars;

fprintf('\n -- ROIs binarised --\n')

%% Adjust segmentation parameters for WF and turn into binary

%%%% Also use ImageJ to find parameters!

minWF = 500;
maxWF = 10000;
contrast(1,1) = 0.01; contrast(1,2) = 0.5;
background(1,1) = 0.01; background(1,2) = 0.5;
binThreshold = 0.7;

close all

% blur the image
I = ICh1;
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

%% save metadata:
metaFile = [name_Ch1 'FOV_' num2str(IM_number) '_Segmentation-Metadata.mat'];
save(metaFile, 'minWF');
save(metaFile, 'maxWF', '-append');
save(metaFile, 'contrast', '-append');
save(metaFile, 'background', '-append');
save(metaFile, 'binThreshold', '-append');

fprintf('\n -- metadata saved --\n')


%% Load localisation data:
close all;

cd(Locpath1);
locs_Ch1 = dlmread(locName1,',',1,0);

% header
file    = fopen(locName1);
line    = fgetl(file);
h       = regexp( line, ',', 'split' );

xCol       = strmatch('x [nm]',h);
yCol       = strmatch('y [nm]',h);
LLCol      = strmatch('loglikelihood',h);

fprintf('\n -- Localisations loaded --\n')


%% Extract rough particle locations from WF - create ROI

box_size = 10;                                                             %diameter of box in pixels

%Find the center of each ROI and transform into an X,Y coordinate 
Center=[]; 
for k=1:length(B)
    boundary    = B{k};
	Center(k,1) = (((max(B{k,1}(:,1))-min(B{k,1}(:,1)))/2)+min(B{k,1}(:,1)))*(pxl);       % Center of the segemented spot in nm
	Center(k,2) = (((max(B{k,1}(:,2))-min(B{k,1}(:,2)))/2)+min(B{k,1}(:,2)))*(pxl);       % Center of the segemented spot in nm
    Center(k,3) = (box_size/2);
end

% Extract the integrated intensity of the WF image for each ROI
intI = []; Particles_WF = [];
for i = 1:length(Center);    
    intI(i,1) = sum(sum(ICh1(int16((Center(i,1)./pxl)-Center(i,3)):int16((Center(i,1)./pxl)+Center(i,3)), int16((Center(i,2)./pxl)-Center(i,3)):int16((Center(i,2)./pxl)+Center(i,3)))));    
    Particles_WF{i,1} = ICh1(int16((Center(i,1)./pxl)-Center(i,3)):int16((Center(i,1)./pxl)+Center(i,3)), int16((Center(i,2)./pxl)-Center(i,3)):int16((Center(i,2)./pxl)+Center(i,3)));
end

% Correction Factor to stretch WF xy dimensions to Loc dimensions
center2 = [];
CFX = (max(locs_Ch1(:,xCol)/pxl))./size(I);
CFY = (max(locs_Ch1(:,yCol)/pxl))./size(I);                                

center2(:,1)    =   Center(:,2)*CFX(:,1);                                  % Center of the segmented spot in nm
center2(:,2)    =   Center(:,1)*CFY(:,1);                                  % Center of the segmented spot in nm
center2(:,3)    =   Center(:,3);                                           % width (=diameter) of cluster
Center          =   center2;

fprintf('\n -- %i ROIs found --\n', length(Center))


% Filter out ROIs with neighbours closer than the box (~500nm)

Distance = []; center2 = [];

count = 1;

for i = 1:length(Center);                                                  % for each ROI
    % Check if radii overlap
    otherClusters = Center(setdiff(1:length(Center),i),:);                 % find all other spots from the list (Center)
    
    for j = 1:length(otherClusters);                                       % for each of these other ROIs
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

fprintf('\n -- Overlapping ROIs filtered: %i ROIs left --\n', length(Center))


% Build box around each Center and copy locs into variable, Cent

Cent1={}; 
count = 1;

for i = 1:length(Center);
    
    vx1 = Center(i,1)+Center(i,3)*pxl;                                   % build box by adding "spot-radius" [note, this was amplified] to center
    vx2 = Center(i,1)-Center(i,3)*pxl;
    vy1 = Center(i,2)+Center(i,3)*pxl;
    vy2 = Center(i,2)-Center(i,3)*pxl;
                                                                           % filter with "spot-boxes":
    target_Ch1 = find(locs_Ch1(:,xCol)>vx2 & locs_Ch1(:,xCol)<vx1 & locs_Ch1(:,yCol)>vy2 & locs_Ch1(:,yCol)<vy1);
    
    Cent1{count,1} = locs_Ch1(target_Ch1,1:end);                            % the localisations
    Cent1{count,2} = length(target_Ch1);                                    % number of localisations
    Cent1{count,3} = intI(i);                                               % WF intensity
    Cent1{count,4} = Particles_WF{i,1};                                     % locations of the box around the particles in the WF         
  
    count = count + 1;   
end

fprintf('\n -- %i Particles selected from localisation dataset --\n',length(Cent1))

% find center of mass of localisations and shift ROIs accordingly
CoM = []; Cent={}; count = 1;

for i = 1:length(Cent1);
    disp(i);
    
    CoM(i,1)    =   median(Cent1{i,1}(:,1));                                % Center of the segmented spot in nm
    CoM(i,2)    =   median(Cent1{i,1}(:,2));
    CoM(i,3)    =   Center(i,3);

% re-build shifted box:
    vx1 = CoM(i,1)+CoM(i,3)*pxl;
    vx2 = CoM(i,1)-CoM(i,3)*pxl;
    vy1 = CoM(i,2)+CoM(i,3)*pxl;
    vy2 = CoM(i,2)-CoM(i,3)*pxl;
                                                                           % filter with "spot-boxes":
    target_Ch1 = find(locs_Ch1(:,xCol)>vx2 & locs_Ch1(:,xCol)<vx1 & locs_Ch1(:,yCol)>vy2 & locs_Ch1(:,yCol)<vy1);
 
    Cent{count,1} = locs_Ch1(target_Ch1,1:end);                           % the localisations
    Cent{count,2} = length(target_Ch1);                                   % number of localisations
    Cent{count,3} = intI(i);                                              % WF intensity
    Cent{count,4} = Particles_WF{i,1};                                    % locations of the box around the particles in the WF         
  
    count = count + 1;  
end

fprintf('\n -- %i Particles selected and from re-aligned WFs. --\n',length(Cent1))

%%
cd(Path);
mkdir('Analysis');
savepath   = [Path 'Analysis'];
cd(savepath)
save(savename,'Cent','-v7.3');

fprintf('\n -- %i Particles saved. --\n',length(Cent1))

%% manually adjust WF-shift to better align: (default = 0, for visualisation purposes only)
xShift = -0.0*pxl;   % >0 moves WF to the right
yShift = -0.0*pxl;   % >0 moves WF down (inverse if want to move locs)

%% Overlay with Widefield image and plot correlations
Cent_wo_duplicates = Cent;

figure('Position',[150 150 400 400],'name',['Extracted Particles from Ch1 on WF']);
imshow(ICh1, [minWF maxWF]); hold on;

for i = 1:size(Cent_wo_duplicates,1);
% backtransform into pxl-coordinates: 1. divide by correction factor + 2. divide through pxl size
    CentO=[];
    CentO(:,1) = Cent_wo_duplicates{i,1}(:,xCol)/CFX(:,1)-xShift;          % divide by correction factor to almost align
    CentO(:,2) = Cent_wo_duplicates{i,1}(:,yCol)/CFY(:,1)-yShift;          % also correct for estimated shift

    scatter(CentO(:,1)/pxl,CentO(:,2)/pxl,1,'magenta');
    hold on;
end

% save overlay:
overlayName = [name_Ch1 'FOV_' num2str(IM_number) '_overlay.fig'];
savefig(overlayName);

% Plot for each Particle the integrate intensity vs. the nbr of locs 
figure('Position',[400 100 300 300],'name','# of Locs vs. WF Intensity');

scatter(cell2mat(Cent_wo_duplicates(:,2)),cell2mat(Cent_wo_duplicates(:,3)),5,'filled','r');hold on;
xlabel('Nbr of locs Ch1');
ylabel('WF int intensity');
box on;
axis square;

fprintf('\n -- Overlay plotted & saved --\n')


end






%% Combine all FOVs into single file
clc;
many = {};

for i = [1,2,3,4];                                                         % can add multiple using 

load([name_Ch1 'FOV_' num2str(i) '_extractedParticles.mat'])

many = vertcat(Cent,many);

clc
fprintf([' Done FOV ' num2str(i)]);

end

save([name_Ch1 'multipleFOVs_all_extractedParticles.mat'],'many','-v7.3');

fprintf('\n - all Finished - ');
