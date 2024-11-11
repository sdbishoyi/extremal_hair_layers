function [Nv, VX, K, EToV] = MeshGen1D_V1(x_regions,K_regions)

% function [Nv, VX, K, EToV] = MeshGen1D(xmin,xmax,K_regions)
% Purpose  : Generate grid, x_regions are boundaries of regions of varying
% number of elements, K_regions is a vector associating number of elements in that
% region

VX=[];
K_total=0;
for ii = 1:max(size(K_regions))
    Nv = K_regions(ii)+1; 
    K_total=K_total+K_regions(ii);
    
    % Generate node coordinates
    VX_temp = (1:Nv);
    for i = 1:Nv
        VX_temp(i) = (x_regions(ii+1)-x_regions(ii))*(i-1)/(Nv-1) + x_regions(ii);
    end
    if(ii~=max(size(K_regions)))
        VX_temp=VX_temp(1:Nv-1);
    end
    VX=[VX VX_temp];
end

% read element to node connectivity
EToV = zeros(K_total, 2);
for k = 1:K_total
    EToV(k,1) = k; EToV(k,2) = k+1;
end

K=K_total;
Nv=K_total+1;

return