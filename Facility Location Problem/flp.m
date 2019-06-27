function [ assignment_vector, Cost ] = flp(problem)
%% FLP - Algorithm to solve the Facility Location Problem 
% Simple version using a greedy strategy, hence does not guarantee the
% optimal solution. However it guarantees a solution in polynomial time.
%
% INPUT:
% - problem: MAT filename with the following variables
%   - Cij: Cost matrix n x m with n individuals and m locations
%   - fij: fixed costs of opening a location m x 1
%
% OUTPUT:
% - assignment_vector: vector with n x 1 size with the location indices
% assigned to each of the n individuals
%
% - Cost: total cost of the solution
%
%
% USAGE EXAMPLE
% [ a_v , Cost ] = flp(problem)
% greedy_assignment = flp(problem)
%
%
%
% Coded by João Araújo, 2019


% Load problem cost matrix Cij and fixed opening costs vector fij
load(problem);

% Initialize variables of interest
u = Inf(size(Cij,1),1);
S = [];
Dij = Cij;
Cost = 0;
initialization = true;

%% Core algorithm
while 1
    
    % Get our minimum Z. If we are not initializing and all Z's > 0, then
    % STOP criterion is reached
    z = sum(Dij,1) + fij';
    [z_min,facility] = min(z);
    
    if z_min > 0
        if initialization
            initialization = false;
        else
            break
        end
    end
    
    % Update facilities set and cost
    S = [S;facility];
    Cost = Cost + z_min;
    
    % Get our u vector from the minimum of this facility costs and previous
    % u vector
    u = min(u,Cij(:,facility));
    
    % Update our Dij matrix as the difference between each column of Cij
    % and u. If this difference is positive then make the element = 0
    Dij = Cij - repmat(u,1,size(fij,1));
    Dij(Dij>0) = 0;
    Dij(:,S) = Inf;

end


% Once we have our best solution, the assigned facilities will be the ones
% with the lowest cost for each individual
S = sort(S);
[~,assignment_vector] = min(Cij(:,S),[],2);
assignment_vector = S(assignment_vector);


end

