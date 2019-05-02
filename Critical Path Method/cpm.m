function solution = cpm(filename,arg2)
% Critical Path Method for activity planning in a project
% 
% INPUT
% filename: name of a file with the following structure:
%           1st col - Name of the activity
%           2nd col - Preceding activities
%           3rd col - time costs of the activities
%
% arg2 (optional): constant string 'schedule' for an extra output of a
%                  project schedule plot
%
%
% OUTPUT
% solution: structure with the following fields
%           - c_path: critical path of the project
%           - M: total time for the project (critical path time)
%           - table: table with activities' earliest times to start (E),
%           latest times to end (L), slack and path from the initial
%           activity to each activity
%
% As an optional feature it can also output the schedule plot. Slacks are
% represented as dashed lines and critical activities as red lines
%
%
%
% USAGE EXAMPLES
% solution = cpm('CPM.dat')
% cpm('CPM.dat','schedule')
% solution = cpm('CPM.dat','schedule')
%
%
%
%
% Coded by João Araújo, 2019
%
%

%% Set-up
% Load raw data from file
fileId = fopen(filename);
Data = textscan(fileId,'%s %s %f');
fclose(fileId);

% Read parameters and build project graph
Graph = buildGraph(Data);

% Load useful variables from graph structure
activities = Graph.struct(:,1);
predecessors = Graph.struct(:,2);
successors = Graph.struct(:,3);
cost = cell2mat(Graph.struct(:,4));

% Initialize all E's, L's and optimal subpaths variables
E = Inf(size(activities,1),1);
L = Inf(size(activities,1),1);
subPaths = cell(size(activities,1),1);

% Convert predecessors and successors to indices for easier manipulation
predecessors = activityToIndex(predecessors,activities);
successors   = activityToIndex(successors,activities);

%% CPM - Earliest times
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
                   subPaths{i} = [subPaths{predecessors{i}(max_idx)},activities{predecessors{i}(max_idx)}];
                    
                end

            end
            
        end

    end
    
end

% Check if we have more than one activitywithout successors. If not, our 
% pseudo-critical path is calculated by adding the last activity. If not, 
% we should do a final iteration of the CPM
if length(final_activities) == 1
    
    M = E(final_activities) + cost(final_activities);
    c_path = [subPaths{final_activities},activities{final_activities}];
    
elseif length(final_activities) > 1
    
   [max_val,max_idx] = max(E(final_activities) + cost(final_activities));

   M = max_val;
   c_path = [subPaths{final_activities(max_idx)},activities{final_activities(max_idx)}];   
    
else
    
    error('Error found in the calculations of Earliest times!')
    
end

%% CPM - Latest times
% Find our latest times knowing that L(activity(i)) =
% min{L(successors(activity(i))) - cost}
while (ismember(Inf,L))
   for i = length(activities):-1:1

        if isempty(successors{i})
            L(i) = E(i) + cost(i);
        else
            
            L_suc = L(successors{i});
            
            if ~any(L_suc == Inf)
                L(i) = min(L_suc - cost(successors{i}));
            end
            
        end
   end
end

% Calculate slack
slack = L - (E+cost);

% Build solution table
sol_table = table(activities,E,L,slack,subPaths,'VariableNames',{'Activities','E','L','slack','Path'});

% Build solution structure
solution.c_path = c_path;     %critical path
solution.M = M;               %critical path time
solution.table = sol_table;   %solution table



%% Schedule plot
% If the user indicates so, generate a plot with the schedule planned for
% the project
if exist('arg2','var') && strcmp(arg2,'schedule')
    
    hold on
    activities_label = cell(size(activities));

    for i = 1:length(activities)
        if slack(i) ~= 0
            plot(E(i):E(i)+cost(i), i * ones(cost(i)+1), 'b-', 'LineWidth',2);
        else
            plot(E(i):E(i)+cost(i), i * ones(cost(i)+1), 'r-', 'LineWidth',2);
        end
        plot(E(i)+cost(i):E(i)+cost(i)+slack(i), i * ones(slack(i)+1), 'b--');
        activities_label{i} = char(activities{i});
    end
    set(gca,'YTick',1:length(activities));
    set(gca,'YTickLabel',activities_label);
    ylabel('Activity');
    xlabel('Time');
    ylim([0,length(activities)+1]);
    xlim([-1,M+1]);
    title('Project Schedule')
    hold off
end




end