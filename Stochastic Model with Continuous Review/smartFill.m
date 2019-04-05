function [prob_struct] = smartFill(prob_struct)
%% Function to calculate problem's missing parameters based on the given data
% Get X and delta from Y and l
if isfield(prob_struct,'Y') && ~isfield(prob_struct,'X')
    Y = prob_struct.Y;

    if strcmp(prob_struct.Y.dist,'uniform') 
        a = prob_struct.Y.a;
        b = prob_struct.Y.b;
        prob_struct.Y.delta = (a+b)/2;
    end
    
    if isfield(prob_struct,'l') && ~isfield(prob_struct,'X')
        l = prob_struct.l;
        if strcmp(prob_struct.Y.dist,'normal')
            delta    = prob_struct.Y.delta;
            sigma_y = prob_struct.Y.sigma;
            X.dist = 'normal';
            X.mu = l * delta; % Because mu = E[X] and mu = l * delta;
            X.sigma = sqrt(l) * sigma_y; % Because sigmax^2 = l*sigmay^2
        elseif strcmp(prob_struct.Y.dist,'uniform')
            X.dist = 'uniform';
            X.a = a * l;
            X.b = b * l; %Because mu = l*delta
            X.mu = (l*(a+b))/2;
        end
        
    % If we do not have 'Y' and 'l' assume l = 1
    else
        prob_struct.l = 1;
        X = Y;
        X.mu = X.delta;
    end
    
    prob_struct.X = X;

    
end

% Get Y from X
if isfield(prob_struct,'X') && ~isfield(prob_struct,'Y')
    
    X = prob_struct.X;
    
    %if we have 'l', we can backtrack Y using previous equations
    if isfield(prob_struct,'l')
        
        if strcmp(prob_struct.Y.dist,'normal')
            mu_x    = prob_struct.X.mu;
            sigma_x = prob_struct.X.sigma;
            Y.dist = 'normal';
            Y.delta = mu_x/l; % Because mu = E[X] and mu = l * delta;
            Y.sigma = sigma_x / sqrt(l); % Because sigmax^2 = l*sigmay^2
        elseif strcmp(prob_struct.Y.dist,'uniform')
            Y.dist = 'uniform';
            Y.a = X.a/l;
            Y.b = X.b/l; %Because mu = l*delta
            Y.delta = (X.a+X.b)/(2*l);
        end
        
    
    % If we do not have 'Y' and 'l' assume l = 1
    else
        
        prob_struct.l = 1;
        Y = X;
        Y.delta = Y.mu;
        
    end
    
    Y.delta = Y.mu;
    prob_struct.Y = Y;
end

% Get l if we have mu and demand per period
if isfield(prob_struct,'X') && isfield(prob_struct,'demand_period') && (~isfield(prob_struct,'l') || prob_struct.l == 1)
    prob_struct.l = prob_struct.X.mu/prob_struct.demand_period;
    l = prob_struct.l;
    if strcmp(prob_struct.Y.dist,'normal')
        mu_x    = prob_struct.X.mu;
        sigma_x = prob_struct.X.sigma;
        Y.dist = 'normal';
        Y.delta = mu_x/l; % Because mu = E[X] and mu = l * delta;
        Y.sigma = sigma_x / sqrt(l); % Because sigmax^2 = l*sigmay^2
    elseif strcmp(prob_struct.Y.dist,'uniform')
        Y.dist = 'uniform';
        Y.a = X.a/l;
        Y.b = X.b/l; %Because mu = l*delta
        Y.delta = (X.a+X.b)/(2*l);
    end
    prob_struct.Y = Y;
end

% Get our maximum quantity to order Q_max if we have delta and maximum time period T_max
if isfield(prob_struct,'T_max')
    prob_struct.Q_max = prob_struct.T_max * prob_struct.Y.delta;
end

% Get our 'r' parameter if we have a defined a minimum safety stock level
if isfield(prob_struct,'safety_stock')
    prob_struct.r = prob_struct.safety_stock + prob_struct.X.mu;
end

end

