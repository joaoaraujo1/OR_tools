function [linkage_dist,agg_sch,C_hist,Delta] = hierarchic(D,method)
%% Hierarchic clustering function with several methods
% Step 1 - Initialize classes array with all the elements belonging to 
% a different class and allocate memory for other variables
C = num2cell(1:numel(D(1,:)));
linkage_dist = nan(length(D)-1,1);
agg_sch = cell(length(D)-1,1);
C_hist = cell(length(D),1);
Delta = nan(size(D));
it_ctr = 1;

% Repeat the algorithm until we only have K classes assigned
while length(C) > 1
    
    % Step 2 - Identify the 2 closest points and create the update indices
    % vector
    min_dist = min(min(D));
    [r,c] = find( D == min_dist );
    r = r(1); c = c(1);
    update_idc = setxor(1:length(D),[r,c]);   
    
    % Step 3 - Update class distances matrix according to the method chosen
    if strcmp(method,'sin') %single linkage (nearest neighbour)
        D(c,update_idc) = min(D([r,c],update_idc));
    elseif strcmp(method,'comp') %complete linkage
        D(c,update_idc) = max(D([r,c],update_idc));
    %elseif strcmp(method,'med') %mean linkage

    %elseif strcmp(method,'ward')%ward method

    end
        
    D(update_idc,c) = D(c,update_idc);
    D = D([1:r-1, r+1:length(D)],[1:r-1, r+1:length(D)]);
    
    %Save our linkage distance
    linkage_dist(it_ctr) = min_dist;
    
    % Update Delta matrix
    Cc = C{c};
    Cr = C{r};

    for i = 1:length(Cc)
        for j = 1:length(Cr)
            Delta(Cc(i),Cr(j)) = linkage_dist(it_ctr);
            Delta(Cr(j),Cc(i)) = Delta(Cc(i),Cr(j));
        end
    end
    
    % Update our agglomeration schedule and history and update classes
    agg_sch{it_ctr} = [C{c},C{r}];
    C_hist{it_ctr} = C;
    it_ctr = it_ctr + 1;
    C(c) = {[C{c},C{r}]};
    C(r) = [];

    
end


end