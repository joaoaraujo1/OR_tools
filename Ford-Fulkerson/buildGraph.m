function Graph = buildGraph(Arcs,vertex_n,Xij,Uij)
%% Build a graph structure with a successors tree and arch properties
% João Araújo, 2019

Graph.struct = cell(vertex_n,2);
Graph.prop = cell(1,2);
idx = 1;

for i = 1:size(Arcs,1)
                    
    % Draw direct arcs where the flow can be increased and indirect arcs
    % where the flow can be decreased
    if Uij(i)-Xij(i) > 0 %direct arc
        temp1 = Graph.struct(Arcs(i,1),1);
        temp2 = Graph.struct(Arcs(i,1),2);
        Graph.struct{Arcs(i,1),1} = [cell2mat(temp1),Arcs(i,2)]; % add node to structure
        Graph.struct{Arcs(i,1),2} = [cell2mat(temp2),1]; % add Direct arc (1)

        Graph.prop{idx,1} = Arcs(i,:); % Arc nodes
        Graph.prop{idx,2} = Uij(i) - Xij(i); % Arc max flow increase
        Graph.prop{idx,3} = 1; % add Direct arc (1)
        
        idx = idx + 1; 
        
    end

    % Draw inverse arc assuming lij = 0
    if Xij(i) > 0
        temp1 = Graph.struct(Arcs(i,2),1);
        temp2 = Graph.struct(Arcs(i,2),2);
        Graph.struct{Arcs(i,2),1} = [cell2mat(temp1),Arcs(i,1)];
        Graph.struct{Arcs(i,2),2} = [cell2mat(temp2),-1]; % add Indirect arch (-1)

        Graph.prop{idx,1} = fliplr(Arcs(i,:));
        Graph.prop{idx,2} = Xij(i);
        Graph.prop{idx,3} = -1; % add Indirect arch (-1)
        
        idx = idx + 1;    

    end
     
end



end