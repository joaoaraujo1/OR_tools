function [l,theta,Paths] = pdm(filename)
%% PDM algorithm for shortest path between a vertex and the network
%
% INPUT: 
% filename: name of the file with the graph properties. The file must have the
% following structure:
% First line (header):
% 1st var: number of vertices
% 2nd var: number of arches
% 3rd var: initial vertex
% Following lines:
% 1st var: arch initial vertex
% 2nd var: arch final vertex
% 3rd var: arch cost
%
%
% OUTPUT:
% l - matrix of history and final cost of shortest path between the initial vertex and vertex i
% theta - matrix of history and final preceding vertex of node i in the shortest path
% Paths - structure for the final solution:
% 1st col: initial and final node
% 2nd col: optimal path
% 3rd col: optimal cost
%
%
%
% USAGE EXAMPLE
%
% [l,theta,Paths] = pdm('pdm2.dat')
%
%
%
% Coded by João Araújo, 2019
%


%% Set-up

% Load data and parameters
Data = dlmread(filename);
Arches = Data(2:end,1:2);
Cij = Data(2:end,end);
vertex_n = Data(1,1);
initial_vertex = Data(1,3);

% Build our graph structure with costs
Graph = buildGraph(Arches,Cij,vertex_n);

% Initialize variables
l = Inf(vertex_n-1,vertex_n); 
l(:,initial_vertex) = 0;
theta = nan(vertex_n-1,vertex_n);

%% Core algorithm

% Initialize K as the set containing the successors of the initial node
K = cell2mat(Graph.struct(initial_vertex,1));
l(1,K) = cell2mat(Graph.struct(initial_vertex,2));
theta(1,K) = initial_vertex;

% Worst case scenario, iterate from 1 to vertex_n - 1
for k = 1:vertex_n-1
    
    % If our k = vertex-1, we have a cycle of negative cost
    if k == vertex_n-1
        error('This graph contains a negative cost cycle!!');
    end
    
    % Initialize new line of l and theta
    l(k+1,:) = l(k,:); 
    theta(k+1,:) = theta(k,:);
    
    % Build successors matrix
    successors_K = [];
    for j = 1:length(K)
        successors_K = [successors_K,cell2mat(Graph.struct(K(j),1))];
    end
    successors_K = unique(successors_K);
    
    % Try to update each successors minimum path
    next_K = [];
    for j = 1:length(successors_K)
        
        % Get our T
        predecessors = cell2mat(Graph.struct(successors_K(j),3));
        predecessors = predecessors( predecessors ~= successors_K(j) );
        T = intersect(predecessors,K);
        
        % Get our cji from every arch T->sucessor
        arches_ji = [T',successors_K(j) * ones(numel(T),1)];
        [~,idx] = ismember(arches_ji,Arches,'rows');
        c_ji = Cij(idx)';
        
        % See if we can get a better minimum path for each successor
        [val,idx] = min([l(k,successors_K(j)),l(k,T) + c_ji]);
        if val ~= l(k,successors_K(j))
            l(k+1,successors_K(j)) = min([l(k,successors_K(j)),l(k,T) + c_ji]);
            next_K = [next_K,arches_ji(idx-1,2)];
            theta(k+1,successors_K(j)) = arches_ji(idx-1,1);
        end
        
    end
    
    % Update K for the next iteration with the vertices that changed
    % minimum path
    K = next_K;
    
    % If there were no updates, terminate
    if sum(l(k+1,:) == l(k,:)) == vertex_n
        break
    end

    
end

%% Calculate shortest paths from thetas

% Truncate l and theta according to our number of iterations
l = l(1:k+1,:);
theta = theta(1:k+1,:);

% Get paths from thetas as indicators of precedence iterating to our
% vertices excluding the initial one
v_set = 1:vertex_n;
v_set = v_set(~ismember(v_set,initial_vertex));
% Initialize Paths structure indicating first and last node, the path
% between them and the total cost
Paths = cell(vertex_n-1,3);
for i = v_set
    
    Paths{i,1} = [initial_vertex,i];
    
    best_path = i;
    pred = NaN; selected = i;
    while(pred ~= initial_vertex)
        pred = theta(end,selected);
        best_path = [best_path,pred];
        selected = pred;
    end
    
    Paths{i,2} = fliplr(best_path);
    Paths{i,3} = l(end,i);
    
end




















end