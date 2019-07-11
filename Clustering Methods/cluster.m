function [C,linkage_dist,agg_sch,Delta,Validation] = cluster(filename,K,arg3,arg4,arg5)
% Cluster - Function for data clustering in K classes using hierarchic methods
%
%
% INPUT
% filename: CSV file with the following structure in each line
%           ID_Number,X_coordinate,Y_coordinate,Z_coordinate...
%
% K: number of classes we want to separate the data into
%
% arg3 - constant string to indicate the method to use
%                   Single Linkage:     'sin'
%                   Complete Linkage:   'comp'
%                   K-means clustering: 'kmeans'
%
% arg4 (optional) - constant string 'plot' to plot de clustered data
%
% arg5 (for kmeans) - maximum number of iterations
%
%
% OUTPUT
% C: Classes containing the ID_numbers of the datapoints.
%
% linkage_dist: history of the linkage distances that were used in each
% clustering step
%
% agg_sch: agglomeration schedule
%
% Delta: final matrix of linkage distances across all elements
%
% Validation: structure with the classification validation measures:
%             Validation.R = cophenetic correlation
%             Validation.S = Stress
%
%
%
% USAGE EXAMPLES
% C = cluster('firstcitiesprime.csv',10,'comp')
% [C,linkage_dist,agg_sch,Delta] = cluster('firstcitiesprime.csv',5,'sin')
% [C,linkage_dist,agg_sch,Delta,Validation] = cluster('firstcitiesprime.csv',10,'sin')
% C = cluster('firstcitiesprime.csv',10,'kmeans',[],100)
%
%
%
% Coded by João Araújo, 2019
%
%



%% Set-up
% Load and read data
Data = csvread(filename);
iD = Data(:,1);
coordinates = Data(:,2:end);

%% Core Algorithm - Hierarchical
if(strcmp(arg3,'sin') || strcmp(arg3,'comp')) %|| strcmp(arg3,'med') || strcmp(arg3,'ward')
    
    % Calculate distance matrix (Euclidean distance)
    D = dist(coordinates,coordinates');
    D(D == 0) = Inf;
    
    [linkage_dist,agg_sch,C_hist,Delta] = hierarchic(D,arg3);
    
    C = C_hist{end-K+1};
    
    %% Validation (Hierarchical): Cophenetic correlation and Stress
    % Initialize our Validation structure
    Validation = struct;

    % Cophenetic correlation
    D(isinf(D)) = 0; Delta(isnan(Delta)) = 0;
    corr = corrcoef(D,Delta);
    R_coph = corr(1,2);

    % Stress measure
    num = 0;
    for i = 1:length(D)-1
        for j = i:length(D)
            num = num + (D(i,j) - Delta(i,j))^2;
        end
    end

    Stress  = num/sum(sum(Delta.^2));

    % Update Validation structure
    Validation.R = R_coph;
    Validation.S = Stress;

%% Core Algorithm - Non-Hierarchical
elseif(strcmp(arg3,'kmeans'))
    
    max_iters = arg5;
    
    % When using K-Means, it is important the initialize the centroids
    % randomly. 
    initial_centroids = kMeansInitCentroids(coordinates, K);
    
    [centroids, C] = runkMeans(coordinates, initial_centroids, max_iters);
    
end



%% Convert class indices in the original iDs and 2D data plot
if(exist('arg4','var') && strcmp(arg4,'plot'))
    subplot(1,2,1)
    plot(coordinates(:,1),coordinates(:,2),'+');
    title('Original Data')
    subplot(1,2,2)
    hold on
end
for k = 1:K
    
    if(exist('arg4','var') && strcmp(arg4,'plot'))
        plot(coordinates(C{k},1),coordinates(C{k},2),'+');
        if strcmp(arg3,'kmeans')
            plot(centroids(k,1),centroids(k,2),'ko');
        end
        C{k} = iD(C{k});
    end
    
end
if(exist('arg4','var') && strcmp(arg4,'plot'))
    title(['Clustered Data with ' num2str(K) ' classes'])
    hold off
end





end