% homemade sign function
function [y]= homeHVSD(x)
[K,N]=size(x);
y = zeros(size(x));

% to left of particle -- hard coded tolerance!
    for i =1:K*(N)
        if(x(i)>=0)
            y(i) = 1;    
        end
    end

end



