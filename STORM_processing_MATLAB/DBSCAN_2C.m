function [DBSCAN_filtered] = DBSCAN_2C(dataDBS, minLength, maxLength, k, Eps, uncertainty);
% Run DBSCAN on each particle 

if isempty(dataDBS) ==1;
    disp('This was empty');
    DBSCAN_filtered        = [];
else
    
    [class,type]           = DBSCAN(dataDBS(:,1:2),k,Eps);
    class2                 = transpose(class);                             % class - vector specifying assignment of the i-th object to certain cluster (m,1)
    type2                  = transpose(type);                              % (core: 1, border: 0, outlier: -1)

    coreBorder             = [];
    coreBorder             = find(type2 >= 0);

    if isempty(coreBorder)==1;
        disp('Borders was empty');
        DBSCAN_filtered    = [];

    else
        subset             = [];
        subset             = dataDBS(coreBorder,1:end);                    % select only core and border points
        subset(:,end+1)    = class2(coreBorder);                           % add column with cluster ID 
        
    % Filter cluster(s):
        DBSCAN_filtered    ={};

        % Determine the ROI-size from raw locs by finding max/min values in each direction:
        roi_size           = [];
        border             = uncertainty;                                  % minimal distance to the ROI-border for any given point of the cluster in [nm] <- set as the max uncertainty
    
        roi_size(1,1)      = max(dataDBS(:,1)) - border;
        roi_size(1,2)      = min(dataDBS(:,1)) + border;
        roi_size(1,3)      = max(dataDBS(:,2)) - border;
        roi_size(1,4)      = min(dataDBS(:,2)) + border;
        
        % Filter    
        for i                  = 1:max(subset(:,end));                     % for all cluster
            vx                 = find(subset(:,end)==i);                   % find the i-th cluster
            selected_Cluster   = [];
            selected_Cluster   = subset(vx,1:end);
        
            % filter clusters that touch the border:
            if selected_Cluster(:,1)<roi_size(1,1) & selected_Cluster(:,1)>roi_size(1,2) & selected_Cluster(:,2)<roi_size(1,3) & selected_Cluster(:,2)>roi_size(1,4);
            % filter clusters by number of localisations:    
            if length(vx)>minLength & length(vx)<maxLength;
            DBSCAN_filtered{size(DBSCAN_filtered,1)+1,1} = selected_Cluster;
            
        %Compute CoM, Rg, Ecc, Length and Width of clusers:
            % Radius of Gyration equals the sum of the variances of x,y(,z) divided by the number of locs 
            DBSCAN_filtered{size(DBSCAN_filtered,1),2}  = sqrt(sum(var(DBSCAN_filtered{size(DBSCAN_filtered,1),1}(:,1:2),1,1)));    

            % Eccentricity corresponds to the covariance of x and y --> sqrt of min/max(Eigenvalues)
            DBSCAN_filtered{size(DBSCAN_filtered,1),3}  = sqrt(max(eig(cov(DBSCAN_filtered{size(DBSCAN_filtered,1),1}(:,1:2))))/min(eig(cov(DBSCAN_filtered{size(DBSCAN_filtered,1),1}(:,1:2)))));
   
            %Length & Width are the max & (perpendicular) Eigenvalues
            DBSCAN_filtered{size(DBSCAN_filtered,1),4}  = sqrt(max(eig(cov(DBSCAN_filtered{size(DBSCAN_filtered,1),1}(:,1:2)))));
            DBSCAN_filtered{size(DBSCAN_filtered,1),5}  = sqrt(min(eig(cov(DBSCAN_filtered{size(DBSCAN_filtered,1),1}(:,1:2)))));
            
            %Centers of Mass (CoM) in [nm]
            DBSCAN_filtered{size(DBSCAN_filtered,1),6} = median(DBSCAN_filtered{size(DBSCAN_filtered,1),1}(:,1));
            DBSCAN_filtered{size(DBSCAN_filtered,1),7} = median(DBSCAN_filtered{size(DBSCAN_filtered,1),1}(:,2));
                    
            % Convex Hull returns the coordinates to describe a convex polyhedron comprising all clustered localisations: DBSCAN_filtered{size(DBSCAN_filtered,1),1}(:,1:2);
            [hull, area] = convhull(DBSCAN_filtered{size(DBSCAN_filtered,1),1}(:,1:2));
            DBSCAN_filtered{size(DBSCAN_filtered,1),8} = hull;
            DBSCAN_filtered{size(DBSCAN_filtered,1),9} = area;
            
            % Alternatively, Convex Hull can be calculated using Delaunay Triangulation:
%             DT = delaunayTriangulation(DBSCAN_filtered{size(DBSCAN_filtered,1),1}(:,1), DBSCAN_filtered{size(DBSCAN_filtered,1),1}(:,2));
%             [Chull, Ahull] = convexHull(DT);
%             DBSCAN_filtered{size(DBSCAN_filtered,1),8} = Chull;
%             DBSCAN_filtered{size(DBSCAN_filtered,1),9} = Ahull;

            end % filter clusters by number of Locs
            end % filter clusters cut by ROI
         end % filter all clusters found
    end % if no clusters found
end % if no input data


