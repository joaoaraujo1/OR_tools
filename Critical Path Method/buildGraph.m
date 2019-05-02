function Graph = buildGraph(Data)
%% Build a graph structure with a successors tree for CPM
% outputs a Graph with the following structure
% col 1 - Activity names
% col 2 - Activity predecessors
% col 3 - Activity successors
% col 4 - Activity time cost
%
% João Araújo, 2019

% Create structures
Graph.struct = cell(size(Data{1},1),4);
successors = cell(size(Data{1},1),1);
predecessors = cell(size(Data{1},1),1);

% Assign parameters
activities = Data{1};
predecessors_data = Data{2};
cost = Data{3};

for i = 1:length(activities)
   
   % Format predecessors and transform them in a successor structure 
    predecessors_it = strsplit(char(predecessors_data(i)),',');
    predecessors_it = predecessors_it(~strcmp(predecessors_it,'-'));
    predecessors{i} = predecessors_it;
    
    for j = 1:length(predecessors_it)
                
        for k = 1:length(activities)
            
            if cell2mat(activities(k)) == cell2mat(predecessors_it(j))
                idx = k;
            end
            
        end
        
        if ~isempty(successors{idx})
            successors{idx} = [successors{idx},activities(i)];
        else
            successors{idx} = activities(i);
        end
        
    end
    
    % Populate the Graph structure with activities
    Graph.struct{i,1} = activities(i);  
    
    % Populate the Graph structure with activity time cost
    Graph.struct{i,4} = cost(i);
    

end

% Populate the graph structure with predecessors
Graph.struct(:,2) = predecessors;

% Populate graph strucutre with successors
Graph.struct(:,3) = successors;




end