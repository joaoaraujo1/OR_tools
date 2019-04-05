function [Q,fi,u,e,F_r] = stochastic(Q,Q0,W,X)
%% Basic iteration of the stochastic model for continuous review

F_r = (1 - Q/W); %FI(u) for normal dist

if strcmp(X.dist,'uniform')
    range = X.b-X.a;
    r = F_r * range;
    fi = (range/2) - r + (r^2)/(2*range);
    u = 0;
elseif strcmp(X.dist,'normal')
    u = norminv(F_r);%fi(r)
    fi = X.sigma * (normpdf(u) - u + u * F_r);
end

Q_new = sqrt(Q0^2 + 2*W*fi);
e = abs(Q-Q_new);
Q = Q_new;


end