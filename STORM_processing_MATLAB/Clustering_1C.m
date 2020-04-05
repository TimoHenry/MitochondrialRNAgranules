%% Filtering of Localisations from segmented ROIs 

%   This script was develeoped by Timo Rey, LEB, IPHYS-EPFL, Switzerland
%                         Final version December 2019

%   Code is based on ParticleFilter-script by Christian Sieben as published in:
%   C. Sieben et al., Nature Methods 15, 777-780 (2018)
%   https://github.com/christian-7
%   This script thus falls under the same licensing as the original code published by C. Sieben [see github].

% Input from roiSegmentation: 
%   1 - locs
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
%   2 - number of locs
%   3 - int Intensity
%   4 - cropped WF ROI

% Workflow: 
%   Load data from single or combined FOVs
%   Filter ROIs based on total number of Locs
%   Filter Locs based on quality
%   Filter Locs based on participation to cluster based on DBSCAN algorithm
%       Computation of size & shape descriptors for each cluster
%   Plot remaining Locs
%   Save as .csv

% Output: 
%   Filtered Localisations as .mat & .csv
%       incl. cluster descriptors
%   Metadata for filtering parameters
%   Image of rendered localisations for each cluster, ordered by size


%% Load the data
clear, clc


acq_Date    =  '20200101';
name_Ch1    =  ['FileName'];

singleFOV   =  0;                                                          % set this to 1 if treating a single FOV:
FOV         =  1;

Path        =  'A:\User\Folder\Analysis\';

outFolder   =  'FolderName';



%%%%%%%%%%%%%%%%%% normally no need to edit below %%%%%%%%%%%%%%%%%%%%%
cd(Path);
batchpath = ['A:\User\Folder'];                                            % location of DBSCAN file

if singleFOV == 1;
    disp('you are treating a single FOV');
    savename   = [name_Ch1 'FOV_' num2str(FOV) '_DBSCAN_filtered.mat'];
    outName    = [name_Ch1 'FOV_' num2str(FOV) '_DBSCAN_filtered.csv'];
    FileName   = [name_Ch1 'FOV_' num2str(FOV) '_extractedParticles.mat'];
    load(FileName);
else
    FOV = 'all';
    disp('you are treating many FOVs');
    savename    =  [name_Ch1 '_DBSCAN_filtered.mat'];
    outName     =  [name_Ch1 '_DBSCAN_filtered.csv'];
    FileName    =  [name_Ch1 '_extractedParticles.mat'];
    load(FileName);
    Cent = many;
end

% re-plot WF vs. # Locs
figure('Position',[400 100 300 300],'name','# of Locs vs. Intensity');
scatter(cell2mat(Cent(:,2)),cell2mat(Cent(:,3)),5,'filled','r');hold on;
xlabel('Nbr of locs');
ylabel('WF integrated intensity');
box on;
axis square;

%% Filter ROIs by the number of locs
MinLocsROI = 2000;    % 2000 for proteins (to get rid of noise) {250 for BrU, 500 for mtDNA}
MaxLocsROI = 15000;   % 25000 for proteins & DNA (to get rid of beads) {15000 for BrU}
count = 0;
Cent_filt1 = {};

for i = 1:length(Cent);    
    if Cent{i,2}>MinLocsROI & Cent{i,2}<MaxLocsROI;
        count = count+1;    

        Cent_filt1{count,1} = Cent{i,1}; 
        Cent_filt1{count,2} = Cent{i,2};
        Cent_filt1{count,3} = Cent{i,3};
        Cent_filt1{count,4} = Cent{i,4};
else end
end

fprintf(['\n' num2str(size(Cent_filt1,1)/size(Cent,1)*100) ' percent left after #Locs filtering with min = ' num2str(MinLocsROI) ' and max = ' num2str(MaxLocsROI) '. \n']);

figure
scatter(cell2mat(Cent_filt1(:,2)),cell2mat(Cent_filt1(:,3)));


%% Quality-filter localisations
close all;

minFrame            = 50;
MinPhotons          = 300;
Maxuncertainty      = 30;
Minsigma            = 100;
Maxsigma            = 400;
startFrame          = 500; % {500 for BrU, 1000 for others}

xCol = 1; yCol = 2; frameCol = 3; uncCol = 4; photonCol = 5; LLCol = 7; sigmaCol = 8; zCol = 12;

for i = 1:length(Cent_filt1)
        
    filter_Ch1          = [];
    filter_Ch1          = find(Cent_filt1{i,1}(:,photonCol) > MinPhotons & Cent_filt1{i,1}(:,uncCol) < Maxuncertainty & Cent_filt1{i,1}(:,frameCol) > minFrame & ... 
                          Cent_filt1{i,1}(:,sigmaCol) > Minsigma & Cent_filt1{i,1}(:,sigmaCol) < Maxsigma & Cent_filt1{i,1}(:,frameCol) > startFrame);

    Cent_filt{i,1}      = Cent_filt1{i,1}(filter_Ch1,1:end);
    Cent_filt{i,4}      = Cent_filt1{i,4};
end

clc
display('  Localizations filtered  ');

%% DBSCAN particle size filter for sweep & batch
count = 1;
tic
% to filter your clusters by # associated localisations:
minSize = 500; % {250 for BrU, 500 for others}
maxSize = 6000; % {4000 for BrU, 6000 for others}
DBSCAN_sweep = {};

cd(batchpath);

for Eps = [12];                                                            % use sequence, f.i. [12,14,16] to test 

    for k = [15];                                                          % can also test multiple
    DBSCAN_filtered = {};

    for m = 1:length(Cent_filt);                                           % can choose subset here
        [DBSCAN_filtered_temp] = DBSCAN_cHull_Edge(Cent_filt{m,1},minSize, maxSize, k, Eps, Maxuncertainty);

        if isempty(DBSCAN_filtered_temp)==1;
        else
            DBSCAN_filtered = vertcat(DBSCAN_filtered,DBSCAN_filtered_temp);
        end
        clear DBSCAN_filtered_temp;
        clc
        X = [' Finished DBSCAN ',num2str(m),' of ',num2str(length(Cent_filt)), ' of Eps = ', num2str(Eps), ' and k = ', num2str(k)];
        disp(X);

        DBSCAN_sweep{count,1} = DBSCAN_filtered;
        DBSCAN_sweep{count,2} = k;
        DBSCAN_sweep{count,3} = Eps;
    end
    count = count +1;    
end
end

cd(Path);

if exist(outFolder);
    cd(outFolder);
else
    mkdir(outFolder);
    cd(outFolder);
end

save(savename,'DBSCAN_sweep','-v7.3');
%save(savename,'DBSCAN_sweep_ConvHull','-v7.3');

fprintf(' -- DBSCAN computed in %f sec -- \n',toc)

%% Save metadata

metaFile = [name_Ch1 '_FOV' num2str(FOV) '_filteringMetadata.mat'];
save(metaFile, 'MinLocsROI');
save(metaFile, 'MaxLocsROI', '-append');
save(metaFile, 'minFrame', '-append');
save(metaFile, 'MinPhotons', '-append');
save(metaFile, 'Maxuncertainty', '-append');
save(metaFile, 'Minsigma', '-append');
save(metaFile, 'Maxsigma', '-append');
save(metaFile, 'startFrame', '-append');
save(metaFile, 'minSize', '-append');
save(metaFile, 'maxSize', '-append');

fprintf('\n -- metadata saved --\n')


%% Save as .csv file
out = [];

for clusterNo = 1:length(DBSCAN_sweep(:,1));
    if isempty(DBSCAN_sweep{clusterNo,1}) == 0;
        out = cell2mat(DBSCAN_sweep{clusterNo,1}(:,2:7));
        out(:,6) = DBSCAN_sweep{clusterNo,2};
        out(:,7) = DBSCAN_sweep{clusterNo,3};
        out(:,8) = cell2mat(DBSCAN_sweep{clusterNo,1}(:,9));
        
        dlmwrite(outName, out, 'delimiter', ',' ,'-append');
    else end
end;

fprintf(' -- DBSCAN saved -- \n')


%% Render the filtered images 

pxlsize = 10.6;                                                            % original size is 106nm -> 10.6 = 10x magnification

current_dir = pwd;
mkdir([name_Ch1 num2str(FOV) '_figures'])
cd([name_Ch1 num2str(FOV) '_figures'])

% set box size:
im_size = []; x = 1; y = 2;
im_size = 80;                                                              % pixel dimensions of rendered (magnified) clusters 
                                                                           % visually inspect clusters do not surpass this

for i = 1:length(DBSCAN_filtered(:,1));

subsetLL = DBSCAN_filtered{i,1};

% Filter the locs:
        heigth = round((max(subsetLL(:,y)) - min(subsetLL(:,y)))/pxlsize);
        width  = round((max(subsetLL(:,x)) - min(subsetLL(:,x)))/pxlsize);
        
        rendered = hist3([subsetLL(:,y),subsetLL(:,x)],[heigth width]);
        empty    = zeros(round(max(max(im_size))*1.5),round(max(max(im_size))*1.5));

        size_DC    = size(empty);
        center_X   = round(size_DC(2)/2);
        center_Y   = round(size_DC(1)/2);

        empty(round(center_Y-heigth/2):(round(center_Y-heigth/2))+heigth-1,round(center_X-width/2):(round(center_X-width/2))+width-1) = rendered;

        name32rendered  = [num2str(round(DBSCAN_filtered{i,2})) '_' num2str(pxlsize) '_ID' num2str(i) '_' acq_Date '.tiff'];% label according to Rg


I32 = [];
I32 = uint32(empty);

t = Tiff(name32rendered,'w');
tagstruct.ImageLength     = size(I32,1);
tagstruct.ImageWidth      = size(I32,2);
tagstruct.Photometric     = Tiff.Photometric.MinIsBlack;
tagstruct.BitsPerSample   = 32;
tagstruct.SamplesPerPixel = 1;
tagstruct.RowsPerStrip    = 16;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.Software        = 'MATLAB';
t.setTag(tagstruct)

t.write(I32);
t.close()

end

cd ..

disp([' -- Images Saved --']);
