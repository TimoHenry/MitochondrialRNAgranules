%% Filtering & Clustering of Localisations from segmented ROIs 

%   This script was develeoped by Timo Rey, LEB, IPHYS-EPFL, Switzerland
%                         Final version March 2020

%   Code is based on ParticleFilter-script by Christian Sieben as published in:
%   C. Sieben et al., Nature Methods 15, 777-780 (2018)
%   https://github.com/christian-7
%   This script thus falls under the same license as the original script by C. Sieben [see github].

% Input from 2C roiSegmentation:
%   1 - locs Channel 1
%       Cent(:,1) [aka 'locs'] contains:
%       1   - x [nm]
%       2   - y [nm]
%       3   - frame
%       4   - uncertainty [nm]
%       5   - intensity [photon]
%       6   - offset [photon]
%       7   - loglikelihood
%       8   - sigma [nm]
%       9   - sigma_y -> = x => = 1
%       10  - dx
%       11  - dy
%   2 - number of locs Channel 1
%   3 - int Intensity Channel 1
%   4 - cropped WF ROI Channel 1
%   5 - 8 analogous to Channel 1, for Channel 2 (d755)

% Workflow: 
%   Load data from single FOVs (to allow overlay with WF)
%   Filter ROIs based on total number of Locs
%   Filter Locs based on quality
%   Filter Locs based on participation to cluster based on DBSCAN algorithm
%       Computation of size & shape descriptors for each cluster
%   Plot remaining Locs
%   Save as .csv

% Output: 
%   Filtered Localisations as .mat & .csv
%       incl. cluster descriptors:
%             Channel 1:
%         1  - Localisations
%             1 - 11 as above
%         2  - Radius of Gyration
%         3  - Eccentricity
%         4  - Length
%         5  - Width
%         6  - CoM x
%         7  - CoM y
%         8  - Indices of points that define ConvexHull polygon
%         9  - ConvexHull area
%         10 - 18 as above for Channel 2:
        
%   Metadata for filtering parameters
%   Image of rendered clusters on WF for each FOV
%   Image of rendered localisations for each cluster, ordered by size


%% Load the data
clear, clc, close all

Path          =  'C:\Users\Public\Documents\03_Architecture\02_STORM_2C\03_processedData\';
FOV           =  12;                                                          

name_Ch     =  ['20200219_TR_COS7_BrU-A647_'];
name_Ch2    =  ['20200219_TR_COS7_D2-d755_'];

outFolder   =  ['20200219_2C_COS7_FOV_' num2str(FOV)];

%%%%%%%%%%%%%%%%%% normally no need to edit below %%%%%%%%%%%%%%%%%%%%%
cd(Path);
% path to DBSCAN-function:
batchpath = ['C:\Users\Public\Documents\03_Architecture\02_STORM_2C\02_ProcessingFiles\MATLAB'];

disp('you are treating a single FOV');
savename   = [name_Ch 'FOV_' num2str(FOV) '_DBSCAN_filtered.mat'];
outName    = [name_Ch 'FOV_' num2str(FOV) '_DBSCAN_filtered.csv'];
FileName   = [name_Ch 'FOV_' num2str(FOV) '_extractedParticles.mat'];
load(FileName);

Cent = Cent5;
% re-plot WF vs. # Locs:
figure('Position',[400 100 300 300],'name','# of Locs vs. Intensity');
scatter(cell2mat(Cent(:,2)),cell2mat(Cent(:,3)),5,'filled','r');hold on;
xlabel('Nbr of locs');
ylabel('WF integrated intensity');
box on;
axis square;

%% Filter ROIs by the number of locs
MinLocs1 = 1000;  % for Channel 1 (to get rid of noise)
MaxLocs1 = 20000; % (to get rid of beads)
MinLocs2 = 1000;  % for Channel 2 (to get rid of noise)
MaxLocs2 = 20000; % (to get rid of beads)

count = 0;
Cent_filt1 = {};

for i = 1:size(Cent,1);    
    if Cent{i,2}>MinLocs1 & Cent{i,2}<MaxLocs1 & Cent{i,6}>MinLocs2 & Cent{i,6}<MaxLocs2;
        count = count+1;    

        Cent_filt1{count,1} = Cent{i,1}; 
        Cent_filt1{count,2} = Cent{i,2};
        Cent_filt1{count,3} = Cent{i,3};
        Cent_filt1{count,4} = Cent{i,4};
        Cent_filt1{count,5} = Cent{i,5}; 
        Cent_filt1{count,6} = Cent{i,6};
        Cent_filt1{count,7} = Cent{i,7};
        Cent_filt1{count,8} = Cent{i,8};

    else end
end

fprintf(['\n' num2str(size(Cent_filt1,1)) ' Foci left after #Locs filtering with min = ' num2str(MinLocs1) ' and max = ' num2str(MaxLocs1) '. \n']);

figure
scatter(cell2mat(Cent_filt1(:,2)),cell2mat(Cent_filt1(:,3)));
%length(Cent_filt1{1,1}(:,1))

%% Quality-filter localisations
% Filter parameters [Ch1 (BrU A647), Ch2 (D2 d750)]

startFrame          = [1, 3500];
MinPhotons          = [300, 300];
Maxuncertainty      = [30, 30];
Minsigma            = [100, 100];
Maxsigma            = [400, 400];

xCol = 1; yCol = 2; frameCol = 3; uncCol = 4; photonCol = 5; LLCol = 7; sigmaCol = 8; zCol = 12;

for i = 1:length(Cent_filt1(:,1))
    
    filter_Ch1          = [];
    filter_Ch1          = find(Cent_filt1{i,1}(:,photonCol) > MinPhotons(1) & Cent_filt1{i,1}(:,uncCol) < Maxuncertainty(1) & Cent_filt1{i,1}(:,frameCol) > startFrame(1) & ... 
                          Cent_filt1{i,1}(:,sigmaCol) > Minsigma(1) & Cent_filt1{i,1}(:,sigmaCol) < Maxsigma(1));

    Cent_filt{i,1}      = Cent_filt1{i,1}(filter_Ch1,1:end);
    Cent_filt{i,4}      = Cent_filt1{i,4};

    % Channel 2:
    filter_Ch2          = [];
    filter_Ch2          = find(Cent_filt1{i,5}(:,photonCol) > MinPhotons(2) & Cent_filt1{i,5}(:,uncCol) < Maxuncertainty(2) & Cent_filt1{i,5}(:,frameCol) > startFrame(2) & ... 
                          Cent_filt1{i,5}(:,sigmaCol) > Minsigma(2) & Cent_filt1{i,5}(:,sigmaCol) < Maxsigma(2));
    disp(i)
    Cent_filt{i,5}      = Cent_filt1{i,5}(filter_Ch2,1:end);
    Cent_filt{i,8}      = Cent_filt1{i,8};  
end

disp('done')
clc
display('  Localizations filtered  ');

%% DBSCAN particle size filter for sweep & batch

count = 1;
tic
% to filter your clusters by # associated localisations:
minSize      = [300,350];
maxSize      = [6000,10000];
DBSCAN_sweep = {};

cd(batchpath);
for Eps = [12]; % can test multiple Eps

    for k = [15]; % can test multiple k
    DBSCAN_filtered = {};

    for m = 1:size(Cent_filt,1); % m == roiID
        
        [DBSCAN_filtered_temp1] = DBSCAN_2C(Cent_filt{m,1}, minSize(1), maxSize(1), k, Eps, Maxuncertainty(1));     %insert k & Eps values here & minLength & maxLength?
        [DBSCAN_filtered_temp2] = DBSCAN_2C(Cent_filt{m,5}, minSize(2), maxSize(2), k, Eps, Maxuncertainty(2));     %insert k & Eps values here & minLength & maxLength?
        if isempty(DBSCAN_filtered_temp1)==1 | isempty(DBSCAN_filtered_temp2)==1 % do not keep clusters if not found one in both channels 
        else
            % check if same dimensions & add 0s if necessary

            if length(DBSCAN_filtered(:,:))>length(DBSCAN_filtered_temp1(:,:))
                % add empty cells:
                deltA = length(DBSCAN_filtered(:,:))-length(DBSCAN_filtered_temp1(:,:));
                where = length(DBSCAN_filtered_temp1(:,:))+deltA;
                DBSCAN_filtered_temp1{size(DBSCAN_filtered_temp1, 1),where} = [];
            else
                deltA = length(DBSCAN_filtered(:,:))-length(DBSCAN_filtered_temp1(:,:));
                where = length(DBSCAN_filtered_temp1(:,:))+deltA;
            end
            % append new clusters from Channel 1:
            DBSCAN_filtered = vertcat(DBSCAN_filtered, DBSCAN_filtered_temp1); 
            
            % find closest clusters in Channel 2 & add:
            CoM_dist = [];
            % Find all distances between clusters:
            for zz = 1:length(DBSCAN_filtered_temp2(:,1)) % for every cluster in Channel 2
                for yy = 1:length(DBSCAN_filtered_temp1(:,1)) % for every cluster in Channel 1
                    % find euclidean distance between centers of mass in x & y
                    CoM_dist(zz,yy) = sqrt((abs(DBSCAN_filtered_temp2{zz,6}-DBSCAN_filtered_temp1{yy,6})+abs(DBSCAN_filtered_temp2{zz,7}-DBSCAN_filtered_temp1{yy,7}))^2);
                end %for every cluster in Ch1
            end %for every cluster in Ch2
            % Find nearest neighbours:
            for zz = 1:length(CoM_dist(:,1))
                [val, yy] = min(CoM_dist(zz,:));
                if val <= min(CoM_dist(:,yy))
                    %it's really the minimum! <- this(zz) Ch2 cluster corresponds to Ch1(yy)
                    [column, row] = find([DBSCAN_filtered{:,6}] == DBSCAN_filtered_temp1{yy, 6}); % find which row in DBSCAN_filtered
                    columns = where-deltA; %length(DBSCAN_filtered(row,:));  %find current length of DBSCAN-filtered

                    % add roiID to column 10 of *temp2:
                    DBSCAN_filtered_temp2{zz,10} = m;
                    
                    for counting = 1:length(DBSCAN_filtered_temp2(zz,:))
                             DBSCAN_filtered(row,columns+counting) = DBSCAN_filtered_temp2(zz,counting);
                    end %counting
                else end %there is a Ch2-cluster closer to that Ch1 cluster
             end %for every cluster (CoM) in Ch2
        end %no clusters found in one of the channels
        X = [' Finished DBSCAN ',num2str(m),' of ',num2str(size(Cent_filt,1)), ' of Eps = ', num2str(Eps), ' and k = ', num2str(k)];
        disp(X);
        
        % Remove clusters from Channel 1 that are not found in Channel 2
        this_row_has_empty = cellfun(@isempty,DBSCAN_filtered); % find empty rows
        DBSCAN_filtered(any(this_row_has_empty,2),:) = []; % delete rows with empty;
        
        DBSCAN_sweep{count,1} = DBSCAN_filtered;
        DBSCAN_sweep{count,2} = k;
        DBSCAN_sweep{count,3} = Eps;
        
    end %end for m
    count = count +1;    
    end %end for k  
end %end for Eps

cd(Path);

if exist(outFolder);
    cd(outFolder);
else
    mkdir(outFolder);
    cd(outFolder);
end

save(savename,'DBSCAN_sweep','-v7.3');

fprintf(' -- DBSCAN computed in %f sec -- \n',toc)

%% Save metadata
metaFile = [name_Ch '_FOV' num2str(FOV) '_filteringMetadata.mat'];
save(metaFile, 'MinLocs1');
save(metaFile, 'MaxLocs1', '-append');
save(metaFile, 'MinLocs2', '-append');
save(metaFile, 'MaxLocs2', '-append');
save(metaFile, 'startFrame', '-append');
save(metaFile, 'MinPhotons', '-append');
save(metaFile, 'Maxuncertainty', '-append');
save(metaFile, 'Minsigma', '-append');
save(metaFile, 'Maxsigma', '-append');
save(metaFile, 'minSize', '-append');
save(metaFile, 'maxSize', '-append');

fprintf('\n -- metadata saved --\n')

%% to plot hull & points, and compute polygon-overlap - save as .csv in the end:
cd(Path); cd(outFolder);
out = []; hullpts = [];
Figures = [name_Ch num2str(FOV) '_figures'];
if exist(Figures);
else
    mkdir(Figures);
end

for clusters = 1:size(DBSCAN_filtered,1)                                   % for every cluster
    disp(num2str(clusters));
    figure('Position',[500 100 500 500],'name',['Convex Hull']);           % make a figure with
    %Channel 2:
    Points2 = DBSCAN_filtered{clusters,10}(:,1:2);                         % spots & convex hull from Ch2 below
    [hu2, area2] = convhull(DBSCAN_filtered{clusters,10}(:,1:2));
    plot(Points2(:,1),Points2(:,2),'g*'); hold on; plot(Points2(hu2,1),Points2(hu2,2),'g-');
    %Channel 1:
    Points1 = DBSCAN_filtered{clusters,1}(:,1:2);                          % Ch1
    [hu, area2] = convhull(DBSCAN_filtered{clusters,1}(:,1:2));
    plot(Points1(:,1), Points1(:,2),'b+'); hold on; plot(Points1(hu,1), Points1(hu,2),'b-');
    
    % save Loc-plots:
    cd(Figures);
    savename_LocPlots = [name_Ch 'FOV_' num2str(FOV) '_' num2str(clusters) '.fig'];
    savefig(savename_LocPlots);

% Then, define an ROI around the cluster & plot Convex Hull. - then calculate the overlap
pNorm1 = []; pNorm2 = [];

    % define ROI-dimensions & crop:
for dimens = [1, 2] % for every dimension:
    %find smallest value of polygon-coordinate for each channel:
    low_limit1 = min(Points1(hu,dimens));
    low_limit2 = min(Points2(hu2,dimens));
    %find smallest between the channels
    if low_limit1 < low_limit2
        low_limit = low_limit1;
    else low_limit = low_limit2;
    end
    %normalise all coordinates by this lower limit:
    pNorm1(:, dimens) = Points1(hu, dimens)-low_limit +40;
    pNorm2(:, dimens) = Points2 (hu2, dimens)-low_limit +40;

    %find largest value of polygon-coordinates for each channel:
    max_limit1 = max(Points1(hu,dimens));
    max_limit2 = max(Points2(hu2,dimens));
    %find largest between the channels
    if max_limit1 > max_limit2
        max_limit = max_limit1;
    else max_limit = max_limit2;
    end
    size_img(1,dimens) = max_limit - low_limit +80;
end
max(size_img);

    % plot the polygons [in nanometer scale] & save these images:
figure('Position',[600 100 500 500],'name',['polygons']);
plot(pNorm1(:,1),pNorm1(:,2),'b-'); hold on; plot(pNorm1(:,1),pNorm1(:,2),'b*');
plot(pNorm2(:,1),pNorm2(:,2),'g-'); plot(pNorm2(:,1),pNorm2(:,2),'g*'); hold off;

    % calculate percentage of overlap
% 1) turn into binary image
img_dim = round(max(size_img));
mask1 = poly2mask(pNorm1(:,1), pNorm1(:,2),img_dim,img_dim);
mask2 = poly2mask(pNorm2(:,1), pNorm2(:,2),img_dim,img_dim);

% Option 1:
% 2) multiply the 2 images -> will leave a binary with only intersect
mask3 = mask1.*mask2;
%imshow(mask3);                                                             % show overlapping region
% 3) calculate areas of binaries & find % of overlap 
area_BrU = nnz(mask1);
area_Overlap = nnz(mask3);

prcnt(clusters,1) = area_Overlap/area_BrU;

fprintf(' -- %f percent of this BrU-cluster is contained withing D2 -- \n',prcnt(clusters,1))

    % calculate euclidean distance between CoMs:
x1 = DBSCAN_filtered{clusters,6};
y1 = DBSCAN_filtered{clusters,7};
x2 = DBSCAN_filtered{clusters,15};
y2 = DBSCAN_filtered{clusters,16};
distance(clusters,1) = pdist2([x1,y1],[x2,y2],'euclidean');

cd(Path); cd(outFolder);
    % Save as .csv:
% Ch1 - define which data:
out(clusters, 1:4)  = cell2mat(DBSCAN_filtered(clusters,2:5));             % cluster descriptors: Rg, Ecc, Length, Width, CoM-x, CoM-y 
out(clusters, 5)  = DBSCAN_filtered{clusters,9};                           % Convex Hull Area
out(clusters, 6)  = prcnt(clusters,1);                                     % Ch1 within Ch2 [overlap] in percent
out(clusters, 7)  = distance(clusters,1);                                  % Distance between CoMs [as fct. of diameter of Ch2?]
% Ch2:
out(clusters, 8)  = DBSCAN_filtered{clusters,11};                          % cluster descriptors: Rg
out(clusters, 9)  = DBSCAN_filtered{clusters,12};                          % Ecc
out(clusters, 10) = DBSCAN_filtered{clusters,13};                          % Length
out(clusters, 11) = DBSCAN_filtered{clusters,14};                          % Width
out(clusters, 12) = DBSCAN_filtered{clusters,18};                          % Convex Hull Area

%%%%%%%%%%%%%%%%%%%%%%
out(clusters, 13) = DBSCAN_filtered{clusters,19};                          % roiID
%%%%%%%%%%%%%%%%%%%%%%

hullpts{clusters,1} = Points1(hu,1);
hullpts{clusters,1}(:,2) = Points1(hu,2);
hullpts{clusters,2} = Points2(hu2,1);
hullpts{clusters,2}(:,2) = Points2(hu2,2);
hullpts{clusters,3} = DBSCAN_filtered{clusters,19};                        % roiID
hullpts{clusters,4} = clusters;                                            % clusterID

fprintf(' -- %f percent of this BrU-cluster is contained withing D2 -- \n',prcnt(clusters,1))

end

% save cluster descriptors
dlmwrite(outName, out, 'delimiter', ',');
fprintf(' -- Clusters saved as .csv -- \n')
% dlmwrite(outName, out, 'delimiter', ',' ,'-append');

%% save hull points
close all
count = 1;
for me = 1:size(hullpts,1)
    for itsLong = 1:size(hullpts{me,1},1)
        %Channel 1
        hullpts2(itsLong, count)   = hullpts{me,1}(itsLong,1);
        hullpts2(itsLong, count+1) = hullpts{me,1}(itsLong,2); 
    end
    for itsLong = 1:size(hullpts{me,2},1)
        hullpts2(itsLong, count+2) = hullpts{me,2}(itsLong,1);
        hullpts2(itsLong, count+3) = hullpts{me,2}(itsLong,2);
    end
    
    hullpts2(:, count+4)     = hullpts{me,3};
    hullpts2(:, count+5)     = hullpts{me,4};

    count = count + 6;
end

hullOut = [name_Ch 'FOV_' num2str(FOV) '_hullPoints.csv'];;
dlmwrite(hullOut, hullpts2, 'delimiter', ',');

fprintf(' -- HullPoints saved -- \n')

%% save Locs:
count = 1; sth = [];

for me = 1:size(DBSCAN_sweep{1,1},1)
    for itsLong = 1:size(DBSCAN_sweep{1,1}{me,1},1)
        %Channel 1
        sth(itsLong, count) = DBSCAN_sweep{1,1}{me,1}(itsLong,1);
        sth(itsLong, count+1) = DBSCAN_sweep{1,1}{me,1}(itsLong,2);
end
    for itsLong = 1:size(DBSCAN_sweep{1,1}{me,10},1)
        % add ConvHull <- indicators for which Loc is part of it.
        %Channel 2
        sth(itsLong, count+2) = DBSCAN_sweep{1,1}{me,10}(itsLong,1);
        sth(itsLong, count+3) = DBSCAN_sweep{1,1}{me,10}(itsLong,2);        
end         
        sth(:, count+4) = me;
    count = count + 5;
end

outName2    = [name_Ch 'FOV_' num2str(FOV) '_DBSCANned_Locs.csv'];
dlmwrite(outName2, sth, 'delimiter', ',' ,'-append');

fprintf(' -- Locs saved as .csv -- \n')

 %% Overlay clusters on WF:
% Find & load WF
pxl = 106;

WF_path         = [Path '../../rawData/' name_Ch2 'WF' num2str(FOV)];
xShift2 = -2.0*pxl;   % >0 moves WF to the right                             % enter estimated shift (from overlay)
yShift2 = -2.0*pxl;   % >0 moves WF down (inverse if want to move locs) 


WF_name        = [name_Ch2 'WF' num2str(FOV) '_MMStack_Pos0.ome.tif'];  
cd(WF_path);
ICh = imread(WF_name);                                                     % Widefield

% plot WF
minWF = 1000; maxWF = 12000;
figure('Position',[150 150 400 400],'name',['Clustered Particles on WF']);
imshow(ICh,[minWF maxWF]); hold on;

LocsToPlot = DBSCAN_filtered;


for i = 1:size(LocsToPlot,1);


% for 2nd color:
Cent2=[];
Cent2(:,1) = LocsToPlot{i,10}(:,xCol)-xShift2;             % divide by correction factor to almost align
Cent2(:,2) = LocsToPlot{i,10}(:,yCol)-yShift2;             % also correct for estimated shift
scatter(Cent2(:,1)/pxl,Cent2(:,2)/pxl,1,'green');

% same for first color:
Cent1=[];
Cent1(:,1) = LocsToPlot{i,1}(:,xCol)-xShift2;             % divide by correction factor to almost align
Cent1(:,2) = LocsToPlot{i,1}(:,yCol)-yShift2;             % also correct for estimated shift
scatter(Cent1(:,1)/pxl,Cent1(:,2)/pxl,1,'blue');

hold on;
end

cd(Path);
cd(outFolder);
savename_OverlayWF = [name_Ch 'FOV_' num2str(FOV) '_overlay_clustersOnWF_inverse.fig'];
savefig(savename_OverlayWF);

fprintf(' -- overlay saved -- \n')


 %% Overlay individual colours on WF:
% Find & load WF
name_Ch1 = name_Ch;

LocsToPlot = DBSCAN_filtered;
pxl = 106;

xShift2 = -2.0*pxl;   % >0 moves WF to the right                             % enter estimated shift (from overlay)
yShift2 = -2.0*pxl;   % >0 moves WF down (inverse if want to move locs) 

% Channel 2
WF_path         = [Path '../../rawData/' name_Ch2 'WF' num2str(FOV)];
WF_name        = [name_Ch2 'WF' num2str(FOV) '_MMStack_Pos0.ome.tif'];  
cd(WF_path);
ICh = imread(WF_name);                                                     % Widefield

    % plot WF
figure('Position',[150 150 400 400],'name',['Clustered D2 on WF']);
imshow(ICh,[minWF maxWF]); hold on;

for i = 1:size(LocsToPlot,1);
Cent2=[];
Cent2(:,1) = LocsToPlot{i,10}(:,xCol)-xShift2;             % divide by correction factor to almost align
Cent2(:,2) = LocsToPlot{i,10}(:,yCol)-yShift2;             % also correct for estimated shift
scatter(Cent2(:,1)/pxl,Cent2(:,2)/pxl,1,'green');

hold on;
end
cd(Path);
cd(outFolder);
savename_OverlayWF = [name_Ch1 'FOV_' num2str(FOV) '_overlay_clustered_D2.fig'];
savefig(savename_OverlayWF);


%% Channel 1
xShift1 = -2.5*pxl;   % >0 moves WF to the right                             % enter estimated shift (from overlay)
yShift1 = -2.5*pxl;   % >0 moves WF down (inverse if want to move locs) 


WF_path         = [Path '../../rawData/' name_Ch1 'WF' num2str(FOV)];
WF_name        = [name_Ch1 'WF' num2str(FOV) '_MMStack_Pos0.ome.tif'];  
cd(WF_path);
ICh = imread(WF_name);                                                     % Widefield

    % plot WF
minWF = 200; maxWF = 2000;
figure('Position',[150 150 400 400],'name',['Clustered BrU on WF']);
imshow(ICh,[minWF maxWF]); hold on;

for i = 1:size(LocsToPlot,1);
% same for first color:
Cent1=[];
Cent1(:,1) = LocsToPlot{i,1}(:,xCol)-xShift1;             % divide by correction factor to almost align
Cent1(:,2) = LocsToPlot{i,1}(:,yCol)-yShift1;             % also correct for estimated shift
scatter(Cent1(:,1)/pxl,Cent1(:,2)/pxl,1,'blue');
hold on;
end


cd(Path);
cd(outFolder);
savename_OverlayWF = [name_Ch2 'FOV_' num2str(FOV) '_overlay_clustered_D2.fig'];
savefig(savename_OverlayWF);

fprintf(' -- overlays saved -- \n')
