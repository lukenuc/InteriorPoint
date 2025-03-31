function [alpha] = LineSearch(v, p_v, tau)

alpha = 1; 
while (v + alpha*p_v < (1-tau)*v)
    alpha = 0.95*alpha;
end

