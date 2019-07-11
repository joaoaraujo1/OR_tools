function idx = findClosestCentroids(X, centroids)
%FINDCLOSESTCENTROIDS computes the centroid memberships for every example
%   idx = FINDCLOSESTCENTROIDS (X, centroids) returns the closest centroids
%   in idx for a dataset X where each row is a single example. idx = m x 1 
%   vector of centroid assignments (i.e. each entry in range [1..K])
%

% Set K
K = size(centroids, 1);

% You need to return the following variables correctly.
idx = zeros(size(X,1), 1);

% ====================== YOUR CODE HERE ======================
% Instructions: Go over every example, find its closest centroid, and store
%               the index inside idx at the appropriate location.
%               Concretely, idx(i) should contain the index of the centroid
%               closest to example i. Hence, it should be a value in the 
%               range 1..K
%
% Note: You can use a for-loop over the examples to compute this.
%

%Tamanho dos dados m
m = size(X,1);

%matriz de todas as distancias dos exemplos aos centroides
distance = zeros(m,K);

%Para cada exemplo calcular a distancia a cada centroide
for i = 1:m
    for k = 1:K
        %defino "Points" como dois pontos: o centroide e o exemplo
        points = [X(i,:); centroids(k,:)];
        %calculo a distancia utilizando a funçao pdist
        distance(i,k) = pdist(points);
    end
    %no final vejo qual e a distancia minima em cada linha e atribuo lhe o
    %centroide consoante o indice onde esta o valor minimo
    [~,idx(i)] = min(distance(i,:));
end
    





% =============================================================

end

