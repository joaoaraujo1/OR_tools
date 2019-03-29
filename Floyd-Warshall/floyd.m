function [lij, pij, lij_hist, pij_hist, nodes_path] = floyd(filename,nodelist)
%% Floyd-Warshall algorithm for longest/shortest paths
%
% INPUT: 
% filename: name of the file with the graph properties. The file must have the
% following structure:
% First line (header) - 4 parameters in order
% 1st var: number of vertices
% 2nd var: number of arches
% 3rd var: path type (0: shortest | 1: longest path)
% 4th var: history (0: verbose | 1: only saves optimal solution)
% Following lines:
% 1st var: arch initial vertex
% 2nd var: arch final vertex
% 3rd var: arch cost
%
% nodelist (optional): specific list of node pairs to automatically calculate the
% longest or shortest paths. dimensions nx2
%
%
% OUTPUT:
% lij - final cost/distances matrix
% pij - final matrix of precedents
% lij_hist - full history of iterations for the distance matrix (verbose
% mode only)
% pij_hist - full history of iterations for the precedents matrix (verbose
% mode only)
% nodespath - structure with the ordered shortest or longest paths and
% respective cost
%
% USAGE EXAMPLE
%
% [l,p] = floyd('FLOYD.DAT')
% [l,p,l_h,p_h] = floyd('FLOYD.DAT')
% [l,p,l_h,p_h,paths] = floyd('FLOYD.DAT',[2,1;2,3;2,4;2,5])
%
%
% Coded by João Araújo, 2019
%


%% Set-up

% Load data and parameters
Data = dlmread(filename);
Edges = Data(2:end,1:end-2);
archCost = Data(2:end,end-1);
vertex_n = Data(1,1);
type = Data(1,3);
noHistory = Data(1,4);
lij_hist = zeros(vertex_n,vertex_n,vertex_n);
pij_hist = zeros(vertex_n,vertex_n,vertex_n);


% Build initial precedence and cost matrices
cij = Inf(vertex_n);
pij = Inf(vertex_n);

% If longest path, replace our Inf with -Inf
if type == 1
    cij = cij * -1;
    pij = pij * -1;
end

% Populate our cost/precedent matrices
for i = 1:vertex_n
    
    pij(i,Edges(Edges(:,1) == i,2)) = i;
    cij(i,Edges(Edges(:,1) == i,2)) = archCost(Edges(:,1) == i);
    
end

lij = cij;

%% Core algorithm

for k = 1:vertex_n % our K is equal to the number of vertices
    
    for i = 1:vertex_n
        for j = 1:vertex_n
            if type == 0 % shortest path
                if lij(i,j) > lij(i,k) + lij(k,j)
                    lij(i,j) = lij(i,k) + lij(k,j);
                    pij(i,j) = pij(k,j);
                end
            elseif type == 1 % longest path
                if lij(i,j) < lij(i,k) + lij(k,j)
                    lij(i,j) = lij(i,k) + lij(k,j);
                    pij(i,j) = pij(k,j);
                end                
            end
        end
        
        % Check for illegal cycles
        if (type == 0 && lij(i,i) < 0) || (type == 1 && lij(i,i) > 0)
            cycle = getcycle(pij,i);
            if type == 0
                errorString = ['The graph contains a negative cost circuit: ' mat2str(cycle)];
            elseif type == 1
                errorString = ['The graph contains a positive cost circuit: ' mat2str(cycle)];
            end
            error(errorString);
        end
        
    end
    
    % update history variables if verbose mode is chosen
    if noHistory == 0
        lij_hist(:,:,k) = lij;
        pij_hist(:,:,k) = pij;
    end
    
end

%% Node pairs path / optional
if exist('nodelist','var')
    nodes_path = {};
    for i = 1:size(nodelist,1)

        nodes_path{i,1} = getpath(pij,nodelist(i,:));
        nodes_path{i,2} = lij(nodelist(i,1),nodelist(i,2));

    end
end













    
    
end