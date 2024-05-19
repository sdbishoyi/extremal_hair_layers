function [y]= homeDirac(n,x)
n=1;
    y = dirac(n,x);
    idx = y == Inf; % find Inf
    y(idx) = 1; 
    idx = y == -Inf; % find Inf
    y(idx) = -1; 

end