function [T,c_path,dist] = pert(filename)
%% Project Evaluation and Review Technique
% PERT technique to estimate how much time a project will take given its
% activities.
%
% INPUT
% filename: name of a file with the following structure:
%           1st col - Name of the activity
%           2nd col - Preceding activities
%           3rd col - a list of values in the format (%f,%f,%f) with the
%           minimum time, mode time and maximum time the activity will take
%           to be completed
%
% OUTPUT:
% T: Time estimate for the whole project with 95% certainty that it won't
% take longer than that value
%
% c_path: Pseudo-critical path for the project
%
% dist: a structure with mean and standard deviation of the normal
% distribution for the possible time the project will take. Includes the
% output of an histogram of this normal distribution
%
% USAGE EXAMPLE:
% [T,c_path,dist] = pert('pert11.dat')
%
%
%
% João Araújo, 2019
%


%% Algorithm
% Load raw data from file
fileId = fopen(filename);
Data = textscan(fileId,'%s %s (%f,%f,%f)');
fclose(fileId);

% Read parameters and build project graph
Graph = buildGraph(Data);

% Search for the pseudo-critical path and get our M and V
[c_path,M] = cpmPert(Graph);
V = sum(cell2mat(Graph.struct(c_path,5)).^2);
activities = Graph.struct(:,1);
c_path = activities(c_path);

% Get our normal distribution N~(M,sqrt(V))
dist.mu = M;
dist.sigma = sqrt(V);

% Get the time where the project will be due with 95% probability, which
% means fi = 1.645 in a normal distribution
T = 1.645 * dist.sigma + dist.mu;

% Simulate a time distribution and output the histogram
n = 10000;
simArr = zeros(1,n);
for i = 1:n
    simArr(i) = normrnd(dist.mu,dist.sigma);
    if simArr(i) < 0
        simArr(i) = 0; % truncated normal
    end
end
hist(simArr,100)



end

