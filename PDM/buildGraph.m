function Graph = buildGraph(Arches,Cij,vertex_n)
%% Build a graph structure with a successors tree (PDM)
% João Araújo, 2019

Graph.struct = cell(vertex_n,3);

for i = 1:size(Arches,1)
    % Successors + cost
    temp1 = Graph.struct(Arches(i,1),1);
    temp2 = Graph.struct(Arches(i,1),2);
    Graph.struct{Arches(i,1),1} = [cell2mat(temp1),Arches(i,2)]; % add node to structure
    Graph.struct{Arches(i,1),2} = [cell2mat(temp2),Cij(i)]; % add cost
    % Predecessors
    temp3 = Graph.struct(Arches(i,2),3);
    Graph.struct{Arches(i,2),3} = [cell2mat(temp3),Arches(i,1)]; % add node to structure
end

end