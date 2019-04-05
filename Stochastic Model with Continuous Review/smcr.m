function [solution_K,solution_N] = smcr(prob_struct,sim_arg,Time)
%% Stochastic Model with Continous Review
% Estimates the parameters of an economic model with stochastic demand and
% simulates a system's behaviour for the indicated period of time
%
% GLOSSARY
% A - fixed cost per delivery
% c - unit cost
% h - housing costs
% p'- cost per unit missing, independent of rupture duration
% Y - demand per unit of time (with known density and distribution)
% r - delivery point
% l - delivery time
% X - demand during delivery time (with known density and distribution)
% a,b - delimiters of the Uniform distribution
% sigma - Normal distribution standard deviation
% mu -  Distribution mean
% T - period duration
% T_max - maximum period duration
% safety_stock - stock stored when the stock reposition is available
% N(Q,r) - quality of service (0-1)
% K(Q,r) - cost of service
% fi(r) = ordinate of the density prob
% FI(r) = Cumulative dist density prob function
%
%
%
% INPUT
% prob_struct: problem structure with the following parameters:
%              - Y,l,A,c,h,p: create a system that minimizes K
%              - Y,l,A,c,h,p,target_serv: create a system that minimizes K
%              with a target quality of service defined
%              - Y,l,A,c,h,p,T_max,safety_stock : Simulate a model that
%              minimizes K and that maximizes N (2 output solutions)
%              - X,p,A,c,h,demand_period: Simulate a model, assuming l = 1
%              and a fixed demand per period of time
%
% sim_arg: optional argument to indicate you want a graphic system
% simulation. If used, should have the string value 'simulate'
%
% Time: For how much time we want the simulation to predict the system's
% behaviour
% 
%
% OUTPUT
% solution_K and solution_N with the following fields:
% r, Q, K, T, c_period[full cost per period], l_serv[level/quality of
% service], I_bar [upper limit for average stock level], safety_stock
%
% USAGE
% solution = smcr(prob_struct)
% [s_k,s_n] = smcr(prob_struct);
% [s_k,s_n] = smcr(prob_struct,'simulate',100);
% 
% 
% 
%
% 

%% Setup

% Try to get other variables from the ones given in the problem
prob_struct = smartFill(prob_struct);

% Load and initialize basic variables
A = prob_struct.A;
h = prob_struct.h;
c = prob_struct.c;
X = prob_struct.X;
Y = prob_struct.Y;
delta = prob_struct.Y.delta;
p = prob_struct.p;
solution_K = struct; solution_N = struct;

% Define stopping criteria
epsilon = 0.001; % minimum improvement to keep iterating
max_it = 10;  % maximum iterations the algorithm will perform

% Define Q0 and W
W = p/h*delta;
Q_0 = sqrt(2*A*delta/h); 


%% Restricted algorithm (restricted T, safety stock or level of service)
if isfield(prob_struct,'T_max') && isfield(prob_struct,'safety_stock') % Time and stock restrictions - get a solution for max N(Q,r) and min K(Q,r)
    
    % Solution for best N(Q,r)
    
    Q = prob_struct.Q_max;
    r = prob_struct.r;
    F_r = (1 - Q/W);
    range = X.b-X.a;
    if strcmp(X.dist,'uniform')
        fi_r = (range/2) - r + (r^2)/(2*range);
    elseif strcmp(X.dist,'normal')
        u = norminv(F_r);%fi(r)
        fi_r = X.sigma * (normpdf(u) - u + u * F_r);
    end
    
    % Value function for costs per unit of time   
    K = A*delta/Q + c*delta + h*(Q/2 + r - X.mu) + (p * delta / Q) * fi_r;
    
    % Costs per period
    I_bar = Q/2 + r - X.mu; %(approximation by excess)
    c_period = (A + c*Q) + (h*(Q/delta)*I_bar) + (p*fi_r);

    % Build solution
    solution_N.r = r;
    solution_N.Q = Q;
    solution_N.K = K;
    solution_N.T = Q/delta;
    solution_N.c_period = c_period;
    solution_N.l_serv = F_r^(delta/Q); %level of service per unit of time
    solution_N.I_bar = I_bar;
    solution_N.safety_stock = r - X.mu;
    
    
    % Solution for best K(Q,r)
    Q = sqrt(Q_0^2 + 2*W*fi_r);

elseif isfield(prob_struct,'target_serv') % Restrictions in the level of service per period
    F_r = prob_struct.target_serv;
    if strcmp(X.dist,'uniform')
        range = X.b-X.a;
        r = F_r * range;
        fi_r = (range/2) - r + (r^2)/(2*range);
    else
        u = norminv(F_r);
        r = u * X.sigma + X.mu;
        fi_r = X.sigma * (normpdf(u) - u + u * F_r);
    end
    %Q = sqrt(Q_0^2 + 2*W*fi_r);
    Q = sqrt(2*A*delta/h)*sqrt((h+p)/p);
    
%% Unrestrained iterative algorithm 
else
    
    Q = Q_0; e = Inf; it = 1;

    while e >= epsilon && it <= max_it

        [Q,fi_r,u,e,F_r] = stochastic(Q,Q_0,W,X);
        it = it + 1;

    end

    if strcmp(X.dist,'uniform')
        r = F_r * (X.b-X.a);
    elseif strcmp(X.dist,'normal')
        r = u * X.sigma + X.mu;
    end

    
end
    
% Value function per unit of time   
K = A*delta/Q + c*delta + h*(Q/2 + r - X.mu) + (p * delta / Q) * fi_r;

% Costs per period
I_bar = Q/2 + r - X.mu; %(approximation by excess)
c_period = (A + c*Q) + (h*(Q/delta)*I_bar) + (p*fi_r);

% Build solution
solution_K.r = r;
solution_K.Q = Q;
solution_K.K = K;
solution_K.T = Q/delta;
solution_K.c_period = c_period;
solution_K.l_serv = F_r^(delta/Q); %level of service per unit of time
solution_K.I_bar = Q/2 + r - X.mu; %(approximation by excess)
solution_K.safety_stock = r - X.mu;


%% System Simulation - run if we have all the necessary info
close all
if exist('sim_arg','var') && strcmp(sim_arg,'simulate') && exist('Time','var');
    T = solution_K.T;
    l = prob_struct.l;
    
    % Setup of Time steps
    is_fraction_T = false;
    if T >= 1
        step = 1; %We can really only simulate for a whole unit of time (as per Y demand)
    else
        step = l;
        disp('WARNING: Cannot accurately perform simulation ( T,l < 1 )')
        is_fraction_T = true;
    end
    
    % Important rule: If E[T] < l, our stock crashes!
    if T < l
        disp('WARNING: T < l, so there will be a stock crash!');
    end
   
    X_arr = step:step:Time;          % time (x-axis)
    Y_arr = zeros(length(X_arr),1);  % stock level (y-axis)
    Y_arr(1) = r;
    r_arr = r*ones(1,length(X_arr)); % r marker


    % Initialize iterator and variable that marks an order is in progress
    i = 1;
    ordering = false;
    while i < length(X_arr)

        % Get our stochastic demand
        if strcmp(Y.dist,'uniform')
            if ~is_fraction_T
                dem = rand * Y.b + Y.a;
            else
                dem = rand * (Y.b*step) + (Y.a*step);
            end
        elseif strcmp(Y.dist,'normal')
            if ~is_fraction_T
                dem = normrnd(Y.delta,Y.sigma);
            else
                dem = normrnd(Y.delta*step,Y.sigma*sqrt(step));
            end
            if dem < 0 % truncate normal
                dem = 0;
            end
        end
        
        % Update stock level
        Y_arr(i+1) = Y_arr(i) - dem; 
        
        % Manage order and update stock if order has arrived
        if Y_arr(i+1) <= r
            if ~ordering
                ordering = true;
                l_count = step;
            elseif ordering && l_count >= prob_struct.l
                ordering = false;
                Y_arr(i+1) = Y_arr(i+1) + Q;
            else
                l_count = l_count + step;
            end

        end
        
        % Increase iterator
        i = i + 1;

    end
    
    %Display plot
    plot(X_arr,Y_arr,'k-');   
    hold on
    xlim([X_arr(1) X_arr(end)]);
    ylim([-Q/10 Q + r]);
    plot(X_arr,r_arr,'b--');
    legend('Stock level','Order point (r)');
    ylabel('Stock units');
    xlabel('Time units');
    hold off
        
    

end









end