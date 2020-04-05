%% DBSCAN function including edge-particle rejection & Convex Hull computation 

%   This script was develeoped by Timo Rey, LEB, IPHYS-EPFL, Switzerland
%                         Final version December 2019

%   Code is based on DBSCAN function by Christian Sieben as published in:
%   C. Sieben et al., Nature Methods 15, 777-780 (2018)
%   https://github.com/christian-7
%   This script thus falls under the same licensing as the original code published by C. Sieben [see github].

% This function allows to:
% - perform DBSCAN clustering on 2D-list of points
% - calculate cluster-descriptors including:
%       Radius of Gyration, Eccentricity, Length & Width, Convex Hull
% - exclude clusters that lie near the edge of the ROI, to avoid arefacts
%       from missing points that lie outside the ROI

%% DBSCAN function:

function [DBSCAN_filtered] = DBSCAN_cHull_Edge(dataDBS,minLength,maxLength,k,Eps,uncertainty)

% Run DBSCAN on each ROI & find clusters:

%k   = 50;                                                                 % minimum number of neighbors within Eps
%Eps = 10;                                                                 % minimum distance between points, nm

if isempty(dataDBS)==1;
    DBSCAN_filtered = [];
    
else
    [class,type]    = DBSCAN(dataDBS(:,1:2),k,Eps);                        % uses parameters specified at input
    class2          = transpose(class);                                    % class - vector specifying assignment of the i-th object to certain cluster (m,1)
    type2           = transpose(type);                                     % (core: 1, border: 0, outlier: -1)

    coreBorder      = [];
    coreBorder      = find(type2 >= 0);

if isempty(coreBorder)==1;
    DBSCAN_filtered = [];
    
else
    subset          = [];
    subset          = dataDBS(coreBorder,1:end);                           % select only core and border points
    subset(:,end+1) = class2(coreBorder);                                  % add column with cluster ID 


% Filter clusters:

DBSCAN_filtered ={};

    % Determine the ROI-size from raw locs by finding max/min values in each direction:
roi_size = [];
border   =  uncertainty; % minimal distance to the ROI-border for any given point of the cluster in [nm] <- set as the max uncertainty
    
roi_size(1,1) = max(dataDBS(:,1)) - border;
roi_size(1,2) = min(dataDBS(:,1)) + border;
roi_size(1,3) = max(dataDBS(:,2)) - border;
roi_size(1,4) = min(dataDBS(:,2)) + border;


for i = 1:max(subset(:,end));                                              % for all clusters
    vx = find(subset(:,end)==i);                                           % find the i-th cluster    

    selected_Cluster = [];
    selected_Cluster = subset(vx,1:end);

    % filter clusters that touch the border:
    if selected_Cluster(:,1)<roi_size(1,1) & selected_Cluster(:,1)>roi_size(1,2) & selected_Cluster(:,2)<roi_size(1,3) & selected_Cluster(:,2)>roi_size(1,4);
    
    % filter clusters by number of localisations:    
    if length(vx)>minLength & length(vx)<maxLength;
        
    DBSCAN_filtered{size(DBSCAN_filtered,1)+1,1} = selected_Cluster;  
    
    
%Compute Rg, Ecc, Length and Width of clusers:
    % Radius of Gyration equals the sum of the variances of x,y(,z) divided by
    % the number of locs
    DBSCAN_filtered{size(DBSCAN_filtered,1),2}  = sqrt(sum(var(DBSCAN_filtered{size(DBSCAN_filtered,1),1}(:,1:2),1,1))); 
    DBSCAN_filtered{size(DBSCAN_filtered,1),7}  = sqrt(sum(var(DBSCAN_filtered{size(DBSCAN_filtered,1),1}(:,1:2),1,1))/2); 

    % Eccentricity 
    % covariance of x and y --> sqrt of min/max(Eigenvalues)
    DBSCAN_filtered{size(DBSCAN_filtered,1),3}  = sqrt(max(eig(cov(DBSCAN_filtered{size(DBSCAN_filtered,1),1}(:,1:2))))/min(eig(cov(DBSCAN_filtered{size(DBSCAN_filtered,1),1}(:,1:2)))));
    
    %Length & Width
    DBSCAN_filtered{size(DBSCAN_filtered,1),4}  = sqrt(max(eig(cov(DBSCAN_filtered{size(DBSCAN_filtered,1),1}(:,1:2)))));
    DBSCAN_filtered{size(DBSCAN_filtered,1),5}  = sqrt(min(eig(cov(DBSCAN_filtered{size(DBSCAN_filtered,1),1}(:,1:2)))));

    % Convex Hull
    % This is the matrix containing all localisations inside this cluster: DBSCAN_filtered{size(DBSCAN_filtered,1),1}(:,1:2);
    [hull, area] = convhull(DBSCAN_filtered{size(DBSCAN_filtered,1),1}(:,1:2));
    DBSCAN_filtered{size(DBSCAN_filtered,1),8} = hull;
    DBSCAN_filtered{size(DBSCAN_filtered,1),9} = area;

    % to plot hull & points:
    % Points = DBSCAN_filtered{size(DBSCAN_filtered,1),1}(:,1:2);
    % [hu, area] = convhull(DBSCAN_filtered{size(DBSCAN_filtered,1),1}(:,1:2));
    % plot(Points(:,1),Points(:,2),'*'); hold on; plot(Points(hu,1),Points(hu,2));
    
    % Convex Hull from Delaunay Triangulation:
    DT = delaunayTriangulation(DBSCAN_filtered{size(DBSCAN_filtered,1),1}(:,1), DBSCAN_filtered{size(DBSCAN_filtered,1),1}(:,2));
    [Chull, Ahull] = convexHull(DT);
    % plot(DT.Points(:,1),DT.Points(:,2),'.','MarkerSize',10); hold on; plot(DT.Points(Chull,1),DT.Points(Chull,2),'r')
    DBSCAN_filtered{size(DBSCAN_filtered,1),10} = Chull;
    DBSCAN_filtered{size(DBSCAN_filtered,1),11} = Ahull;    
    
    else end
    else end
    
end

end
end

end


