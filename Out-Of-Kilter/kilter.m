function [Xij,reduced_cost,omegas,total_cost] = kilter(filename,solution)
%% Out-Of-Kilter algorithm for minimum cost flow 
%
% INPUT: 
% filename: name of the file with the graph properties. The file must have the
% following structure:
% First line (header) - 4 parameters in order
% 1st var: number of vertices
% 2nd var: number of arches
% Following lines:
% 1st var: arch initial vertex (i)
% 2nd var: arch final vertex (j)
% 3rd var: arch minimum capacity (lij)
% 4th var: arch maximum capacity (uij)
% 5th var: arch cost (cij)
%
% solution (optional): structure with an initial solution. Should have the
% fields solution.Xij and solution.omegas
% Possible initial solution for the graph 'dados3.ou1':
% solution.Xij = [3 2 3 5 5]';
% solution.omegas = [2 0 -3 2]';
%
%
% OUTPUT:
% Xij: final array of flows for each arch (Primal vector)
% reduced_cost: final reduced cost Cij - omega(i) + omega(j)
% omegas: Dual vector
% total_cost: total cost of the solution 
%
% USAGE EXAMPLE
%
% [Xij,reduced_cost,omegas,total_cost] = kilter('dados2.ou1')
% [Xij,reduced_cost,omegas,total_cost] = kilter('dados3.ou1',solution)
%
%
% Coded by João Araújo, 2019
%

%% Set-up

% Load data and parameters
Data = dlmread(filename);
vertex_n = Data(1,1);
vertices = 1:vertex_n;
arch_n = Data(1,2);
Arches = Data(2:end,1:2);
Lij = Data(2:end,3); % minimum flow
Uij = Data(2:end,4); % maximum flow
Cij = Data(2:end,5); % Cost


% Initialize out-of-kilter table data
if exist('solution','var')
    Xij = solution.Xij;
    omegas = solution.omegas;
    reduced_cost = Cij - omegas(Arches(:,1)) + omegas(Arches(:,2));
else
    Xij = zeros(arch_n,1);
    omegas = zeros(vertex_n,1);
    reduced_cost = Cij; 
end
kilter_state = false(arch_n,1);
total_cost = 0;


%% Core algorithm
while(1)
        
    % Update kilter states
    kilter_state = updateStates(kilter_state,Xij,omegas,Lij,Uij,Cij,Arches);
    
    % build the list of arches we will try to get in-kilter
    ook_list = Arches(kilter_state == false,:);
    
    % We have out-of-kilter arches. Try primal solution
    if any( kilter_state == false) 
        
        % Get auxiliary graph with out-of-kilter arches and use their upper
        % limit as cost (adm)
        Graph = buildGraph(Arches,vertex_n,Xij,reduced_cost,Uij,Lij);
        
        % Check if we have a cycle containing an out-of-kilter arch
        [breakthrough,X0] = scanCycle(Graph,ook_list);

        % If no breakthrough was achieved, procede to dual solution
        if ~breakthrough
            X0 = unique(cell2mat(X0));
            X0_bar = setxor(X0,vertices);
            
            % Get indices for S1 and S2
            S1_idx = find(ismember(Arches(:,1),X0) & ismember(Arches(:,2),X0_bar) & reduced_cost > 0 & Xij <= Uij);
            S2_idx = find(ismember(Arches(:,2),X0) & ismember(Arches(:,1),X0_bar) & reduced_cost < 0 & Xij >= Lij);
            S_idx = sort([S1_idx;S2_idx]);
            
            if isempty(S_idx)
                error('S1 and S2 are empty sets!!')
            end

            % Calculate theta          
            theta = min(abs(reduced_cost(S_idx)));
            
            % If theta is Inf, there is no admissible solution
            if theta == Inf
                error('This problem does not have a solution!!')
            end
            
            % Update omegas
            omegas(X0) = omegas(X0) + theta;
            % Get new dual solution
            reduced_cost = Cij - omegas(Arches(:,1)) + omegas(Arches(:,2));
        
        % If we achieved a breakthrough update flows (primal solution)
        else
            
            % Get our arches, arch type (direct or indirect) and calculate
            % delta as the maximum capacity of the lowest capacity arch
            C = cell2mat(X0(:,1));
            delta = min(cell2mat(X0(:,2)));
            C_type = cell2mat(X0(:,3));
            
            % Update our xij for our direct arches (C_type == 1)
            [~,idx] = ismember(C(C_type == 1,:),Arches,'rows');
            Xij(idx) = Xij(idx) + delta;
            
            % Update our xij for our indirect arches (C_type == -1)
            [~,idx] = ismember(fliplr(C(C_type == -1,:)),Arches,'rows');
            Xij(idx) = Xij(idx) - delta;
            
        end

    % All arches are in-kilter. Get total cost and terminate
    else
        total_cost = sum(Cij.*Xij);
        break        
    end

end

    
end