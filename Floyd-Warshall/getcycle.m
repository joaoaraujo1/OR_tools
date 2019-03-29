function cycle = getcycle(pij,idx)
%% Get a known cycle from a precendents matrix table

ptr = idx;
cycle = idx;
while(1)
    
    node = pij(idx,ptr);
    cycle = [cycle,node];
    ptr = node;
    
    if node == idx
        break
    elseif length(cycle) == length(pij)
        error('No cycle starting/ending with the specified index can be found');
    end

end

cycle = fliplr(cycle);

end