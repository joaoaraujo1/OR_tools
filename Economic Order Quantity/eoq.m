function [solution] = eoq(prob_struct)
%% Economic Order Quantity
% Function that applies the EOQ to a problem. Features:
% - Basic EOQ
% - EOQ with quantity discounts
% - EOQ with stock rupture (w/instant stock reposition, R = Inf)
% - EOQ with stock rupture and quantity discounts (numeric/iterative method)
%
% INPUTS
% prob_struct: problem structure with the following fields:
%              - demand/time (D)
%              - fixed cost/order (A)
%              - housing costs/unit/time (h)
%              - stock break cost/unit/time (p) [only for stock rupture variant of eoq]
%              - fixed cost (c) [basic/stock rupture variant] or array of costs [only for quantity discounts variant]
%              - Minimum quantity Q(i) to order to get the price c(i) [only for quantity discounts variant]
%              - (scale) for numeric Q* estimation (only for stock rupture + quantity discounts variant)
%
%
% OUTPUTS
% solution: structure with the following fields:
%           - Q_star: Optimal quantity to order
%           - T_star: Optimal time interval to order
%           - K_star: K value for Q_star
%           // specific for variants with stock rupture
%           - S_max_star: Optimal maximum lacking stock
%           - I_max_star: Optimal maximum stock housed
%           - S_bar: Average lacking stock/period
%           - I_bar: Average housed stock /period
%           
%
%
% Coded by João Araújo, 2019
%

% Load variables
D = prob_struct.D;
A = prob_struct.A;
h = prob_struct.h;
cQ = prob_struct.c;


% Check if we want a quantity discount version of the eoq
quantity_discount = false;
if isfield(prob_struct,'Q')
    Q = prob_struct.Q;
    quantity_discount = true;
else
    Q = 1;
end

% Check if we want a stock rupture version of the eoq
rupture_allowed = false;
if isfield(prob_struct,'p')
    if prob_struct.p ~= Inf
        p = prob_struct.p;
        rupture_allowed = true;
    end
end

% Calculate Q~
if rupture_allowed
    Q_tilda = sqrt(2*A*D/h)*(sqrt((h + p)/p));
else
    Q_tilda = sqrt(2*A*D/h);
end

% Basic EOQ or just one of its variants
if ~rupture_allowed || ~quantity_discount

    % Get indices array from where we calculate our possible Q*
    indices = find(Q >= Q_tilda);

    % If no index >= Q~ calculate Q* for our closest index / only entry
    if isempty(indices)
        cQ_best = cQ(length(cQ));

        if rupture_allowed
            K_star = cQ_best*D + sqrt(2*A*D*h) * sqrt(p/(h+p));
        else
            K_star = A*D/Q_tilda + cQ_best*D + h*Q_tilda/2;
        end

        Q_star = Q_tilda;

    % Otherwise, calculate all K's from the Q~ index onwards [quantity discounts variant] 
    else

        K_star = Inf;

        for i = indices

            K_temp = A*D/Q(i) + cQ(i)*D + h*Q(i)/2;

            % If our calculated K is smaller than our best K, overwrite our best K
            if K_temp < K_star
                K_star = K_temp;
                Q_star = Q(i);
            end

        end

    end

    % Get our T* from Q*
    T_star = Q_star / D;

    % Build solution for basic eoq and eoq with discounts in quantity
    solution.Q_star = Q_star;
    solution.K_star = K_star;
    solution.T_star = T_star;

    % Build solution for the stock rupture variant of the model
    if rupture_allowed
        solution.S_max_star = sqrt(2*A*D*h/(p*(h+p)));
        solution.I_max_star = Q_star - solution.S_max_star;
        solution.S_bar = solution.S_max_star^2 / (2*Q_star);
        solution.I_bar = solution.I_max_star^2 / (2*Q_star);
    end

% Model with a mix of the 2 EOQ variants (numeric/graphic resolution)    
else
    scale = prob_struct.scale;
    
    %Check if our last Q in the quantity discounts is larger than the
    %minimum of a simple EOQ value function
    if Q(length(Q)) < Q_tilda
        QQ = 1:scale:Q_tilda+1; 
    else
        QQ = 1:scale:Q(length(Q));
    end
    
    nq = length(QQ);
    KK = Inf(nq,1);
    
    %Iterate across all our possible Q's and calculate K, knowing the
    %relationship between Q* and Smax*
    for i = 1:length(QQ)
        Q_temp = QQ(i);
        S_temp = Q_temp * h/(h+p); % Mandatory relationship between Q* and Smax* -> Smax*/Q* = h/(h+P)
        idx = find(Q <= Q_temp,1,'last');
        KK(i) = A*D/Q_temp + cQ(idx)*D + h/2*Q_temp - h*S_temp + ((h+p)/2)*(S_temp^2/Q_temp);
    end
    [~,idx] = min(KK);
    Q_star = QQ(idx); K_star = KK(idx);
    
    %Build approximate solution for the mixed variants model
    solution.Q_star = Q_star;
    solution.K_star = K_star;
    solution.T_star = Q_star / D;
    solution.S_max_star = Q_star*h/(h+p);
    solution.I_max_star = Q_star - solution.S_max_star;
    solution.S_bar = solution.S_max_star^2 / (2*Q_star);
    solution.I_bar = solution.I_max_star^2 / (2*Q_star);
    
    % Plot cost function of K
    subplot(1,2,1)
    semilogy(QQ,KK,'b');
    hold on
    plot(QQ(idx),KK(idx),'r.','MarkerSize',10);
    xlabel('Q'); ylabel('K(Q,Smax(Q))');
    title('K value across Q where Smax(Q) = Q * h/(h+p)')
    legend('Value(K)','Optimal Solution');
    
    % EXTRA: Create 3D grid to get the solution space of K with every
    % possible Q and Smax near the calculated optimum
    SS = solution.S_max_star/1.5:solution.S_max_star*1.5; ns = length(SS);
    QQ = solution.Q_star/1.5:solution.Q_star*1.5; nq = length(QQ);
    KK = Inf(nq,ns);
    for i = 1:length(QQ)
        for j = 1:length(SS)
            Q_temp = QQ(i); S_temp = SS(j); 
            idx = find(Q <= Q_temp,1,'last');
            KK(i,j) = A*D/Q_temp + cQ(idx)*D + h/2*Q_temp - h*S_temp + ((h+p)/2)*(S_temp^2/Q_temp);
        end
    end
    subplot(1,2,2)
    surf(SS,QQ,KK); 
    xlabel('Smax'); ylabel('Q'); zlabel('K(Q,Smax)')
    title('Solution space of K across Q and Smax')
    
end


end