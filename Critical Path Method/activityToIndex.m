function activity_tree = activityToIndex(activity_tree,activities)
%% Activity to Index
% Function to convert predecessor and successor structures in more workable
% numeric indices
%
% João Araújo, 2019

for i = 1:length(activity_tree)
    
    conv_idx = [];
    
    for j = 1:length(activity_tree{i})
        
        for k = 1:length(activities)
            
            idx = find(cell2mat(activity_tree{i}(j)) == cell2mat(activities{k}),1);
            
            if ~isempty(idx)
                conv_idx = [conv_idx,k];
                idx = [];
                break
            end
            
        end

    end
    
    if ~isempty(conv_idx)
        activity_tree{i} = conv_idx;
    end
    
end


end