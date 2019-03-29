function [has_cycle, X0, fringe] = scanCycle(Graph,ook_list)
%% Check if there is a cycle in our auxiliary graph
% Chooses the first out-of-kilter arch of the list and explores the search
% space using a greedy strategy. Search is controlled by a fringe variable
% and a closed set of already explored arches.
%
% João Araújo,2019
%

% Initialize variables
has_cycle = false;
arches_aux = cell2mat(Graph.prop(:,1));
maxflows = cell2mat(Graph.prop(:,2));
arches_type = cell2mat(Graph.prop(:,3));

% Choose the arch to test for cycle and update out-of-kilter list
critical_arch = ook_list(1,:);
%critical_arch = ook_list(randi(size(ook_list,1)),:); % random out-of-kilter node selection. slower but can lead to different dual/primal solutions (with same minimum cost). Try it!
[~,idx] = ismember(critical_arch,arches_aux,'rows');

% If our critical_arch is saturated (meaning it is now an indirect arch), reverse it 
if idx == 0
    critical_arch = fliplr(critical_arch);
    [~,idx] = ismember(critical_arch,arches_aux,'rows');
end

goal = critical_arch(1); % if we reach this node we have a cycle
X0{1,1} = critical_arch;
X0{1,3} = arches_type(idx);
closed_set = [critical_arch(1),critical_arch(1)]; % add arch to the closed set to avoid double arch cycle

% Initialize search tree structure
fringe = [];
candidates = cell2mat(Graph.struct(critical_arch(1),1));
candidates = candidates(~ismember(candidates,closed_set));

for i = 1:size(candidates,2)
    [~,idx] = ismember([critical_arch(1),candidates(i)],arches_aux,'rows');
    
    entry_path = [critical_arch(1), candidates(i)];
    entry_path_flow = maxflows(idx);
    entry_maxflow = maxflows(idx);
    entry_type = arches_type(idx);
    
    fringe = [fringe;{entry_path},{entry_path_flow},entry_maxflow,entry_maxflow,{entry_type}];
end

% Core search algorithm
while((length(closed_set) < length(arches_aux)) && ~isempty(fringe) && ~has_cycle)
    
    % Greedy search: Choose the highest flow arch
    [~,arch_idx] = max(cell2mat(fringe(:,3)));
    critical_path = cell2mat(fringe(arch_idx,1));
    critical_node = critical_path(end);
    candidates = cell2mat(Graph.struct(critical_node,1));
    candidates = candidates(candidates ~= critical_path(end-1));
    
    if ~isempty(candidates) % we have successor nodes to explore
        candidates = [ones(size(candidates,2),1)*critical_node,candidates'];
        candidates = candidates(~ismember(candidates,closed_set,'rows'),:);
        for i = 1:size(candidates,1)
            [~,idx] = ismember(candidates(i,:),arches_aux,'rows');
            
            entry_path = [critical_path,candidates(i,2)];
            
            past_archflows = cell2mat(fringe(arch_idx,2));
            entry_path_flow = [past_archflows, maxflows(idx)];
            
            entry_maxflow = sum(entry_path_flow);
            
            past_type = cell2mat(fringe(arch_idx,5));
            
            cand_type = arches_type(idx);
            entry_type = [past_type,cand_type];
            
            [~,indices] = unique(entry_path);
            
            if candidates(i,2) == goal % we have a cycle
                has_cycle = true;
                fringe = [{entry_path},{entry_path_flow},entry_maxflow,entry_maxflow,{entry_type}];
                break; 
            elseif candidates(i,2) ~= goal && numel(indices) ~= numel(entry_path) % we have a subcycle
                indices = setdiff(1:size(entry_path, 2), indices);
                [~,indices2] = unique(fliplr(entry_path));
                indices2 = setdiff(1:size(entry_path, 2), indices2);
                indices2 = length(entry_path) - (indices2-1);
                entry_path = entry_path(indices2:indices);
                entry_path_flow = entry_path_flow(indices2:indices-1);
                entry_maxflow = entry_maxflow - maxflows(idx);
                entry_type = entry_type(indices2:indices-1);
                has_cycle = true;
                fringe = [{entry_path},{entry_path_flow},entry_maxflow,entry_maxflow,{entry_type}];
                break; 
            end
            fringe = [fringe;{entry_path},{entry_path_flow},entry_maxflow,entry_maxflow,{entry_type}];
        end
        
        if ~has_cycle
            fringe(arch_idx,:) = []; %delete fringe from previous state;
        end
        
    else % deadend
        fringe(arch_idx,3) = {-1}; 
    end
    
    closed_set = [closed_set;closed_set(end),critical_node];

end

% Retrieve the nodes we can reach or cycle as a set of arches
if ~isempty(fringe)
    it = 1;
    X0 = {[0,0],[],[]};
    for k = 1:size(fringe,1)
        
        nodes = cell2mat(fringe(k,1));

        for i = 1:size(nodes,2)-1
            node_to_add = [nodes(i),nodes(i+1)];
            if ~ismember(node_to_add,cell2mat(X0(:,1)),'rows')
                [~,arch_idx] = ismember(node_to_add,arches_aux,'rows');
                X0{it,1} = [nodes(i),nodes(i+1)];
                X0{it,2} = abs(maxflows(arch_idx));
                X0{it,3} = arches_type(arch_idx);
                it = it + 1;
            end
        end
    end
    
    % If we did not find a cycle, exclude the initial node from X0
    if ~has_cycle
        X0 = unique(cell2mat(X0(:,1)));
        X0 = {X0(X0~=goal)};
    end
else % if we can only reach a node our starting point, pass that node in X0
    X0 = {critical_arch(2)};
end





end