function [c_path,M] = cpmPert(Graph)
%% CPM method for usage in the PERT algorithm
% Uses a algorithm to find the Earliest time (E) for every
% activity and outputs the pseudo-critical path as well as the total
% estimated cost of that path
%
%
% João Araújo, 2019

%% Set-up
% Load useful variables from graph structure
activities = Graph.struct(:,1);
predecessors = Graph.struct(:,2);
successors = Graph.struct(:,3);
cost = cell2mat(Graph.struct(:,4));

% Initialize all E's and subpaths
E = Inf(size(activities,1),1);
subPaths = cell(size(activities,1),1);

% Convert activity names to indices for easier manipulation
for i = 1:length(predecessors)
    
    conv_idx = [];
    
    for j = 1:length(predecessors{i})
        
        for k = 1:length(activities)
            
            idx = find(cell2mat(predecessors{i}(j)) == cell2mat(activities{k}),1);
            
            if ~isempty(idx)
                conv_idx = [conv_idx,k];
                idx = [];
                break
            end
            
        end

    end
    
    if ~isempty(conv_idx)
        predecessors{i} = conv_idx;
    end
    
end

%% Core algorithm
final_activities = [];

% Find all the Earliest dates for each activity
while (ismember(Inf,E))

    for i = 1:length(activities)
        
        if isempty(successors{i}) && ~ismember(i,final_activities)
           final_activities = [final_activities;i];
        end
            
          
        
        % Update all the E's that were not yet calculated
        if E(i) == Inf
            
            % If there are no predecessors, the Earliest start is 0
            if isempty(predecessors{i})

                E(i) = 0;

            else
                
                % If we have the E's of all the predecessors we can
                % calculate the E for the present activity
                if sum(E(predecessors{i}) ~= Inf) == length(predecessors{i})
                    
                   [max_val,max_idx] = max(E(predecessors{i}) + cost(predecessors{i}));
                   
                   E(i) = max_val;
                   subPaths{i} = [subPaths{predecessors{i}(max_idx)},predecessors{i}(max_idx)];
                    
                end

            end
            
        end

    end
    
end

% Finally, to get the full cost, check if we have more than one activity
% without successors. If not, our pseudo-critical path is calculated by
% adding the last activity. If not, we should do a final iteration of the
% CPM
if length(final_activities) == 1
    
    M = E(final_activities) + cost(final_activities);
    c_path = [subPaths{final_activities},final_activities];
    
elseif length(final_activities) > 1
    
   [max_val,max_idx] = max(E(final_activities) + cost(final_activities));

   M = max_val;
   c_path = [subPaths{final_activities(max_idx)},final_activities(max_idx)];   
    
else
    
    error('Error found in the CPM algorithm')
    
end





end