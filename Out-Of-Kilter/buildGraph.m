function Graph = buildGraph(Arches,vertex_n,Xij,reduced_cost,Uij,Lij)
%% Build a graph structure with a successors tree and arch properties
% João Araújo, 2019

Graph.struct = cell(vertex_n,2);
Graph.prop = cell(1,3);
idx = 1;

for i = 1:size(Arches,1)
                    
    % Draw direct arch and its flow according to the in-Kilter table
    if ((reduced_cost(i) < 0 || reduced_cost(i) == 0) && Xij(i) < Uij(i)) || Xij(i) < Lij(i) %S1 excluding 0 cost arches
        temp1 = Graph.struct(Arches(i,1),1);
        temp2 = Graph.struct(Arches(i,1),2);
        Graph.struct{Arches(i,1),1} = [cell2mat(temp1),Arches(i,2)]; % add node to structure
        Graph.struct{Arches(i,1),2} = [cell2mat(temp2),1]; % add Direct arch (1)

        Graph.prop{idx,1} = Arches(i,:); % Arch nodes
        
        if Xij(i) < Lij(i) && reduced_cost(i) >= 0
            Graph.prop{idx,2} = Lij(i) - Xij(i); % Arch maxflow
        else
            Graph.prop{idx,2} = Uij(i) - Xij(i); % Arch maxflow
        end
        
        Graph.prop{idx,3} = 1; % add Direct arch (1)
        
        idx = idx + 1; 
    end

    % Draw inverse arch according to the in-Kilter table
    if ((reduced_cost(i) > 0 || reduced_cost(i) == 0) && Xij(i) > Lij(i)) || Xij(i) > Uij(i) %S2 excluding 0 cost arches
        temp1 = Graph.struct(Arches(i,2),1);
        temp2 = Graph.struct(Arches(i,2),2);
        Graph.struct{Arches(i,2),1} = [cell2mat(temp1),Arches(i,1)];
        Graph.struct{Arches(i,2),2} = [cell2mat(temp2),-1]; % add Indirect arch (-1)

        Graph.prop{idx,1} = fliplr(Arches(i,:));
        
        if Xij(i) > Uij(i) && reduced_cost <= 0
            Graph.prop{idx,2} = Uij(i) - Xij(i);
        else
            Graph.prop{idx,2} = Xij(i) - Lij(i);
        end
        
        Graph.prop{idx,3} = -1; % add Indirect arch (-1)
        
        idx = idx + 1;    

    end
     
end



end