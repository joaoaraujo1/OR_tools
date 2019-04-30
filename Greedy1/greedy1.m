function [ items_list, max_value ] = greedy1( problem_table )
% Greedy1 heuristic to solve the Knapshack Problem (maximization)
%
% INPUT:
% problem_table - a table with the linear problem with a structure similar
% to the one in the example (tableau)
%
% OUTPUT
% items_list: a list of the best items to include in the knapshack
% according to the heuristic
% max_value: value of the solution
%
% USAGE EXAMPLES
% ---> load('mat_with_problem_table') previously
%
% [items_list, max_value] = greedy1(problem_table)
% items_list = greedy1(problem_table)
% 
%
%
% Coded by João Araújo, 2019
%

%% Set-up
% Read problem and get our parameters
items = problem_table.Properties.VariableNames;
values = table2array(problem_table({'value'},1:end-1));
weights = table2array(problem_table({'weights'},1:end-1));
capacity = table2array(problem_table({'weights'},end));

% Sort our variables by value/weight ratio
[~,item_idx] = sort(values./weights,'descend');


%% Core algorithm
% While the knapshack is not at full capacity try to insert the items in 
% descending order of value/weight
cur_weight = 0;
items_n = 0;
items_list = {};
max_value = 0;
for i = 1:size(weights,2)
    
    % If we have space, insert item
    if cur_weight + weights(item_idx(i)) <= capacity
        
        items_list(items_n+1) = items(item_idx(i));
        items_n = items_n + 1;
        max_value = max_value + values(item_idx(i));
        cur_weight = cur_weight + weights(item_idx(i));
        
        % We reached our maximum capcity. End algorithm
        if cur_weight == capacity
            break
        end

    end
    
end





end

