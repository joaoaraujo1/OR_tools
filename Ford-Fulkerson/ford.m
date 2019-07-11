function [ max_flow, flow_vec, X0, X0_bar, min_cut ] = ford( filename )
%% FORDFULKERSON Algorithm to calculate the maximum feasable flow on a network
%  
% INPUT:
% - filename: name of the file with the following structure
%             1st line - number of vertices TAB number of arcs TAB initial
%             vertex TAB final vertex
%             following lines - arc initial vertex TAB arc final vertex TAB arc maximum capacity (uij)
%
% OUTPUT
% - maxflow: the maximum flow allowable in the network
% - flow_vec: a structure with the arcs, flow and maximum capacity
% that result in the best solution
% - X0 and X0_bar: X0 are the nodes one can reach from the initial node
% after the final solution is applied. X0_bar is the complementary set
% - min_cut: minimum cut to disconnect the network
%
% USAGE EXAMPLES
% [ max_flow, flow_vec, X0, X0_bar, min_cut ] = ford('FLUMAX.dat')
% max_flow = ford('FFSLIDES.dat')
%
%
% Coded by João Araújo
%

%% Set-up

% Load data and parameters
Data = dlmread(filename);
vertex_n = Data(1,1);
vertices = 1:vertex_n;
arc_n = Data(1,2);
initial_vertex = Data(1,3);
final_vertex = Data(1,4);
verbose = Data(1,5);
Arcs = Data(2:end,1:2);
nnz_lij = false;
max_flow = 0;

if max(Data(2:end,3)) > 0 && max(Data(2:end,4)) > 0 
    nnz_lij = true; % non-zero lower limit for some arcs. Apply special steps. *NOT IMPLEMENTED YET*
    Lij = Data(2:end,3); % minimum flow
    Uij = Data(2:end,4); % maximum flow
    error('Current implementation does not support lij>0')
else
    Uij = Data(2:end,3); % maximum flow
end

Xij = zeros(arc_n,1);

%% Core algorithm
while(1)
    
    % Build auxiliary graph
    Graph = buildGraph(Arcs,vertex_n,Xij,Uij);
    
    % Search for a path v(initial->final)
    [found_path, X0] = scanPath(Graph,initial_vertex,final_vertex);
    
    % If a path was found, update flows Xij
    if found_path
        
        % Get our arcs, arc type (direct or indirect) and calculate
        % delta as the capacity of the lowest capacity arch
        C = cell2mat(X0(:,1));
        delta = min(cell2mat(X0(:,2)));
        C_type = cell2mat(X0(:,3));

        % Update our xij for our direct arcs (C_type == 1)
        [~,idx] = ismember(C(C_type == 1,:),Arcs,'rows');
        Xij(idx) = Xij(idx) + delta;

        % Update our xij for our indirect arcs (C_type == -1)
        [~,idx] = ismember(fliplr(C(C_type == -1,:)),Arcs,'rows');
        Xij(idx) = Xij(idx) - delta;
        
        max_flow = max_flow + delta;
        
    % If no path was found, we have reached our max solution. Calculate X0
    % and X0_bar as well as the minimum cut
    else
        X0_bar = setxor(X0,vertices);
        
        % Get indices for X0->X0_bar and X0_bar->X0
        S1_idx = find(ismember(Arcs(:,1),X0) & ismember(Arcs(:,2),X0_bar));
        S2_idx = find(ismember(Arcs(:,2),X0) & ismember(Arcs(:,1),X0_bar));
        S_idx = sort([S1_idx;S2_idx]);
        min_cut = Arcs(S_idx,:);
        flow_vec = cell(size(Xij,1),3);
        
        % Build flow vector with Uij and arc info
        for i = 1:size(Xij,1)
            flow_vec(i,1) = {Arcs(i,:)};
            flow_vec(i,2) = {Xij(i)};
            flow_vec(i,3) = {Uij(i)};
        end
        
        break
         
    end
    

    
    
end



end

