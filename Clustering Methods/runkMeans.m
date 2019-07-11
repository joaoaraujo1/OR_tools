function [centroids, C] = runkMeans(X, initial_centroids, max_iters)
%RUNKMEANS runs the K-Means algorithm on data matrix X, where each row of X
%is a single example
%   [centroids, idx] = RUNKMEANS(X, initial_centroids, max_iters, ...
%   plot_progress) runs the K-Means algorithm on data matrix X, where each 
%   row of X is a single example. It uses initial_centroids used as the
%   initial centroids. max_iters specifies the total number of interactions 
%   of K-Means to execute.r unkMeans returns 
%   centroids, a Kxn matrix of the computed centroids and idx, a m x 1 
%   vector of centroid assignments (i.e. each entry in range [1..K])
%

% Initialize values
[m n] = size(X);
K = size(initial_centroids, 1);
centroids = initial_centroids;
idx = zeros(m, 1);

% Run K-Means
for i=1:max_iters
    
    % For each example in X, assign it to the closest centroid
    idx = findClosestCentroids(X, centroids);
    
    % Given the memberships, compute new centroids
    centroids = computeCentroids(X, idx, K);
end

% Convert idx in C class structure
C = cell(K,1);
for k = 1:K
    C{k} = find(idx == k);
end


end

