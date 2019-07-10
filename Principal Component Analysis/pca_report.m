function [ descriptives, lambdas, Aij, Rij, Lij, Yij, explained, var_weights, relevant_vars, relevant_pcs] = pca_report( Dataset, variance_threshold, alpha_var )
%% PCA_REPORT - A complete multivariate exploratory analysis with PCA
%   
% INPUTS:
% - Dataset: name of the CSV without headers containing the raw dataset
% - variance_threshold: allowed proportion of change between maximum and 
% minimum variance variables before we standardize the dataset
% - alpha_var: alpha for dimension reduction criterion 1
%
%
% Coded by João Araújo, 2019

standardize = false;
Xij = csvread(Dataset);

%% Descriptives

descriptives.mean = nanmean(Xij); % mean
descriptives.std = nanstd(Xij); % standard deviation
descriptives.var = nanvar(Xij); % variance
descriptives.cov = nancov(Xij); % covariance
descriptives.corr = corrcoef(Xij); % pearson correlation
descriptives.gen_var = det(descriptives.cov); % generalized variance
descriptives.tot_var = sum(descriptives.var); % total variance


%% Decide if dataset should be standardized

if min(descriptives.var) < (1-variance_threshold) * max(descriptives.var)
    standardize = true;
    Xij = (Xij - repmat(descriptives.mean,size(Xij,1),1))./repmat(descriptives.std,size(Xij,1),1);   
    fprintf('\nVariances too dissimilar. Using standardized dataset\n\n');
else
    fprintf('\nSimilar variances. Using raw dataset\n\n');
end

%% PCA: eigenvectors/eigenvalues, correlations and loadings

% Get lambdas and full PC matrix (Aij) using SVD
[Aij,Diaglambda] = svd(cov(Xij));
lambdas = diag(Diaglambda);

% Get Xij-Aij correlations (Rij)
if standardize
    Rij = Aij.*repmat(sqrt(lambdas'),size(Aij,1),1);
else
    Rij = Aij.*repmat(sqrt(lambdas'),size(Aij,1),1)./repmat(descriptives.std,size(Aij,1),1);
end

% Get our loadings matrix (Lij)
if standardize
    Lij = Rij;
else
    Lij = Rij.*repmat(descriptives.std,size(Aij,1),1);
end

% Get our scores (Yij)
Yij = Xij*Aij;

%% PCA: metrics
% Proportion of original variance explained by each PC
explained.global = lambdas./sum(lambdas);

% Proportion of variance of each variable explained by each PC
explained.variables = Lij.^2;

% Relative weight of each variable for each principal component
var_weights = Aij.^2; % Since sum(aij^2) = 1

% Decide the relevant variables for each PC based on 2 criteria
% Criterium 1: |Rij| >= .5
% Criterium 2: |Rij| >= sqrt(rij^2/p)
relevant_vars.c1 = {}; relevant_vars.c2 = {};
for j = 1:size(Xij,2)
    arr_c1 = find(abs(Rij(:,j)) >= .5);
    arr_c2 = find(abs(Rij(:,j)) >= sqrt(sum(Rij(:,j).^2)/size(Xij,1)));
    relevant_vars.c1 = [relevant_vars.c1;{arr_c1}];
    relevant_vars.c2 = [relevant_vars.c2;{arr_c2}];
end

% Decide the relevant PCs for the dataset based on 2 criteria
% Criterium 1: Predefined minimum variance explained (based on the input
% alpha_var)
% Criterium 2: Variance explained needs to be larger than the mean of the
% variance explained by all PCs
alpha = alpha_var;
lambdas_cum = zeros(size(lambdas));
for j = 1:size(Xij,2)
    lambdas_cum(j) = sum(lambdas(1:j)/size(Xij,2));
end

mean_exp = sum(lambdas)/size(Xij,2); % If dataset is standardized, this value is 1, so we are using Kaiser rule
relevant_pcs.c1 = find(lambdas_cum >= alpha); relevant_pcs.c1 = 1:min(relevant_pcs.c1);
relevant_pcs.c2 = find(lambdas(lambdas >= mean_exp));










end

