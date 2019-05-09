function C = knn(filename,K,arg3)
% KNN - Nearest Neighbours algorithm for data clustering in K classes
%
%
% INPUT
% filename: CSV file with the following structure in each line
%           ID_Number,X_coordinate,Y_coordinate,Z_coordinate...
%
% K: number of classes we want to separate the data into
%
% arg3 (optional): constant string 'plot' to output a 2D plot of the
%                  original and clustered data
%
%
% OUTPUT
% C: Classes containing the ID_numbers of the datapoints.
%
% If used with arg3 will also output the data plots
%
%
% USAGE EXAMPLES
% C = knn('firstcitiesprime.csv',10)
% C = knn('firstcitiesprime.csv',5,'plot')
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

% Calculate distance matrix (Euclidean distance)
D = mandist(coordinates,coordinates');
D(D == 0) = Inf;

%% Core Algorithm

% Step 1 - Initialize classes array with all the elements belonging to a different
C = num2cell(1:numel(D(1,:)));

% Repeat the algorithm until we only have K classes assigned
while length(C) > K
    
    % Step 2 - Assign the 2 closest points to the same class and delete the
    % empty class
    min_dist = min(min(D));
    [r,c] = find( D == min_dist );
    r = r(1); c = c(1);
    update_idc = setxor(1:length(D),[r,c]);
    C(c) = {[C{c},C{r}]};
    C(r) = [];
    
    % Step 3 - Update class distances matrix
    D(c,update_idc) = min(D([r,c],update_idc));
    D(update_idc,c) = D(c,update_idc);
    D = D([1:r-1, r+1:length(D)],[1:r-1, r+1:length(D)]);
    
    
end

%% Convert class indices in the original iDs and 2D data plot
if(exist('arg3','var') && strcmp(arg3,'plot'))
    subplot(1,2,1)
    plot(coordinates(:,1),coordinates(:,2),'+');
    title('Original Data')
    subplot(1,2,2)
    hold on
end
for k = 1:K
    
    if(exist('arg3','var') && strcmp(arg3,'plot'))
        plot(coordinates(C{k},1),coordinates(C{k},2),'+');
    end
    C{k} = iD(C{k});
    
end
title('Clustered Data')
hold off





end