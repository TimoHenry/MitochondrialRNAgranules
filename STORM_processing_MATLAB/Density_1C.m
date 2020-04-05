%% Density estimations of clusters

% INPUT: .m-files with cluster descriptors
% OUTPUT: .csv-files with cluster descriptors, including number of locs
% WORKFLOW:
    % *DBSCAN_filtered.m files to load cluster-data
    % extract number of locs from first column by ConvexHull Area [column 9]
    % save new .csv-file with number of locs.

%% Load the data
clear, clc, close all

Date          = '20200101';
name          = 'FileName_multipleFOVs_all_DBSCAN_filtered.mat';

Probe         = '1hBrU';
Cells         = 'COS7';
outName       = [Cells '_' Probe '_' Date '.csv'];

Path_part1    =  'A:\User\Folder\';
Path_part2    = [Path_part1 Cells '_' Probe '_' Date '\Analysis\filtered_MRGs\'];

cd(Path_part2);
load(name);
%% save as .csv
out = [];
for clusterNo = 1:length(DBSCAN_sweep(:,1));
    if isempty(DBSCAN_sweep{clusterNo,1}) == 0;
        out = cell2mat(DBSCAN_sweep{clusterNo,1}(:,2:7));
        out(:,6) = DBSCAN_sweep{clusterNo,2};
        out(:,7) = DBSCAN_sweep{clusterNo,3};
        out(:,8) = cell2mat(DBSCAN_sweep{clusterNo,1}(:,9));
        for i = 1:size(DBSCAN_sweep{1,1},1)
            out(i,9) = size(DBSCAN_sweep{clusterNo,1}{i,1},1);
        end
        dlmwrite(outName, out, 'delimiter', ',' ,'-append');
    else end
end;
fprintf(' -- Clusters saved incl. #Locs -- \n')

cd(Path_part1);