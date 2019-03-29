function path = getpath(pij,node_tuple)
%% Get a known path from a precendents matrix table

startNode = node_tuple(1);
endNode = node_tuple(2);
ptr = endNode;
path = endNode;
while(1)
    
    node = pij(startNode,ptr);
    path = [path,node];
    ptr = node;
    
    if node == startNode
        break
    elseif length(path) == length(pij) || node == Inf
        error('No path can be found between specified nodes');
    end

end

path = fliplr(path);

end