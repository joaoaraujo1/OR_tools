function [found_path, X0] = scanPath(Graph,initial_vertex,final_vertex)
%% Check for maximum v(initial_node->final_node) in our auxiliary graph
% Greedy searcher
% fringe structure (per line):
% 1) Entry path - nodes in the current search path
% 2) Entry path flow - max flows of each arc in the current search path
% 3) and 4) maximum flow possible. 3) gets -1 value when a deadend is reached
% 5) entry_type - array of 1's (direct) and -1's (indirect) according to the arc type
%
% João Araújo,2019
%

% Initialize variables
found_path = false;
arcs_aux = cell2mat(Graph.prop(:,1));
maxflows = cell2mat(Graph.prop(:,2));
arcs_type = cell2mat(Graph.prop(:,3));
first_iteration = true;
fringe = [];
closed_set = [NaN,NaN];

% Core search algorithm
while((length(closed_set) <= length(arcs_aux)) && ~found_path)
    
    % If iteration = 1 start search with initial node
    % If iteration > 1 Greedy search: Choose the max flow arc
    if first_iteration
        critical_node = initial_vertex;
        candidates = cell2mat(Graph.struct(critical_node,1));
        critical_path = critical_node; arc_idx = 1;
    else
        [~,arc_idx] = max(cell2mat(fringe(:,3)));
        critical_path = cell2mat(fringe(arc_idx,1));
        critical_node = critical_path(end);
        candidates = cell2mat(Graph.struct(critical_node,1));
        candidates = candidates(candidates ~= critical_path(end-1));
    end
    
    if ~isempty(candidates) % we have successor nodes to explore
        candidates = [ones(size(candidates,2),1)*critical_node,candidates'];
        candidates = candidates(~ismember(candidates,closed_set,'rows'),:);
        for i = 1:size(candidates,1)
            [~,idx] = ismember(candidates(i,:),arcs_aux,'rows');
            
            entry_path = [critical_path,candidates(i,2)];
            
            if ~first_iteration
                past_archflows = cell2mat(fringe(arc_idx,2));
            else
                past_archflows = [];
            end
            entry_path_flow = [past_archflows, maxflows(idx)];
            entry_maxflow = min(entry_path_flow);
            
            if ~ first_iteration
                past_type = cell2mat(fringe(arc_idx,5));
            else
                past_type = [];
            end
            
            cand_type = arcs_type(idx);
            entry_type = [past_type,cand_type];
                        
            if candidates(i,2) == final_vertex % we have a path
                found_path = true;
                fringe = [{entry_path},{entry_path_flow},entry_maxflow,entry_maxflow,{entry_type}];
                break; 
            end
            fringe = [fringe;{entry_path},{entry_path_flow},entry_maxflow,entry_maxflow,{entry_type}];
        end
        
        if ~found_path && ~first_iteration
            fringe(arc_idx,:) = []; %delete fringe from previous state;
        end
        
    else % deadend
        if ~isempty(fringe)
            fringe(arc_idx,3) = {-1}; 
        else % isolated first node
            break
        end
    end
    
    closed_set = [closed_set;closed_set(end),critical_node];
    if first_iteration
        first_iteration = false;
        closed_set = closed_set(2:end,:);
    end

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
                [~,arc_idx] = ismember(node_to_add,arcs_aux,'rows');
                X0{it,1} = [nodes(i),nodes(i+1)];
                X0{it,2} = abs(maxflows(arc_idx));
                X0{it,3} = arcs_type(arc_idx);
                it = it + 1;
            end
        end
    end
    
    % If we did not find a path, exclude the initial node from X0
    if ~found_path
        X0 = unique(cell2mat(X0(:,1)));
    end
else % if we can only reach a node - our starting point - pass that node in X0
    X0 = initial_vertex;
end





end