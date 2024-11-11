% A purely outgoing ell = 2 initial data for
%
%        [\partial_t^2 - \partial_r + ell(ell+1)/r^2]psi = 0
%
% is given by 
%
%         psi(t,r) = f'' + (3/r)f' + (3/r^2) f
%
% for some profile function f(u) of retarded time u = t-r
%
%
% For nargout == 3:
% This routine generates outgoing data (t=0) for a sine-gaussian profile 
% 
% f(u) = sin[f_0(u-x0)] exp(-c(u-x0)^2) 
% 
% 
% For nargout == 1:
% This routine generates the recorded time-series psi(t-r1,r2)
%
%  When using hyperboloidal layers, call this function as follows: 
%
%    Psi_exact_inf = OutgoingEll2_ID_SineGaussian(rend,end),f0,c,x0,simulation_time_grid - x(end,end));
%
% where..
%   x is the rho grid 
%   r is the usual radial coordiante 
%   simulation_time_grid are the times to evaluate the solution at (tau)

function [Psi,Psi_r,Psi_t,Psi_tt] = OutgoingEll2_ID_SineGaussian(x,f0,c,x0,u)

           %%% nargout = 4 -- compute initial data %%%
                 
% Input: computational grid x
%        sine wave frequence f0
%        Gaussian width c
%        Gaussian center x0
%        retarded time u = t-x = tau - rho (in hyperboloidal coordinates) 
%
%
%  This can be used in different ways
%
% which can be used to get initial data
%
%  [Psi,Psi_r,Psi_t]=OutgoingEll2_ID_SineGaussian(x,f0,c,x0,0)
%
% the solution at some later time (say t=30)
%
%  [Psi,Psi_r,Psi_t]=OutgoingEll2_ID_SineGaussian(x,f0,c,x0,30)
%
%
% or the solution at some location on the grid (say r=30)
%
%  [Psi,Psi_r,Psi_t]=OutgoingEll2_ID_SineGaussian(30,f0,c,x0,0:.1:100)
%
%
%
%
%
%

if(nargout == 1)
    
    %%% Compute derivaties of f on the shifted grid X %%%
    %u = time - x;
    X = u - x0;
    
    %f = sin(f0*X).*exp(-c*X.^2);
    %f1 = (f0*cos(f0*X) - 2*c*X.*sin(f0*X)).*exp(-c*X.^2);
    %f2 = ( (4*(c^2).*X.^2 - 2*c - f0^2).*sin(f0*X) - 4*c*f0.*X.*cos(f0*X) ).*exp(-c*X.^2);
    %f3 = ( 2*c.*X.*(-4*(c^2).*X.^2 + 6*c + 3*f0^2).*sin(f0*X) - f0*(-12*(c^2).*(X.^2)+6*c+f0^2).*cos(f0*X)).*exp(-c*X.^2);
    
    %%% Assemble outgoing initial data %%%
    % psi(t,r) = f'' + (3/r)f' + (3/r^2) f and psi_t and psi_r from chain rule
    f2 = ( (4*(c^2).*X.^2 - 2*c - f0^2).*sin(f0*X) - 4*c*f0.*X.*cos(f0*X) ).*exp(-c*X.^2);

    Psi = zeros(max(size(X)),1); 
    Psi(:) = f2(:);
    
    
end

if(nargout == 4)
    
    %%% Compute derivaties of f on the shifted grid X %%%
    %u = time - x;
    X = u - x0;
    
    f = sin(f0*X).*exp(-c*X.^2);
    f1 = (f0*cos(f0*X) - 2*c*X.*sin(f0*X)).*exp(-c*X.^2);
    f2 = ( (4*(c^2).*X.^2 - 2*c - f0^2).*sin(f0*X) - 4*c*f0.*X.*cos(f0*X) ).*exp(-c*X.^2);
    f3 = ( 2*c.*X.*(-4*(c^2).*X.^2 + 6*c + 3*f0^2).*sin(f0*X) ...
         - f0*(-12*(c^2).*(X.^2) ...
                + 6*c ...
                + f0^2).*cos(f0*X)).*exp(-c*X.^2);
    f4 = exp(-c*X.^2).*( ...
        (16*c^4*(X.^4) - 48*c^3 * (X.^2) + c^2 * (12 - 24*f0^2 *(X.^2)) + 12*c*f0^2 + f0^4).*sin(f0*X) ...
        - 8*c*f0.*X.*(4*c^2.*(X.^2) - 6*c - f0^2).*cos(f0*X));

    Psi_tt = f4 + (3./x).*f3 + (3./x.^2).*f2;
    
    %%% Assemble outgoing initial data %%%
    % psi(t,r) = f'' + (3/r)f' + (3/r^2) f and psi_t and psi_r from chain rule
    
    Psi    = f2 + (3./x).*f1 + (3./x.^2).*f;
    Psi_t  = f3 + (3./x).*f2 + (3./x.^2).*f1;
    Psi_r  = -f3 - (3./x.^2).*f1 - (3./x).*f2 - (3./x.^2).*f1 - (6./x.^3).*f;

    %%% Diagnostic: check violation of compact support %%%
    XTemp = x(end-1):.01:2*x(end-1);
    X = -XTemp - x0;
    
    f = sin(f0*X).*exp(-c*X.^2);
    f1 = (f0*cos(f0*X) - 2*c*X.*sin(f0*X)).*exp(-c*X.^2);
    f2 = ( (4*(c^2).*X.^2 - 2*c - f0^2).*sin(f0*X) - 4*c*f0.*X.*cos(f0*X) ).*exp(-c*X.^2);
    f3 = ( 2*c.*X.*(-4*(c^2).*X.^2 + 6*c + 3*f0^2).*sin(f0*X) - f0*(-12*(c^2).*(X.^2)+6*c+f0^2).*cos(f0*X)).*exp(-c*X.^2);


    PsiTemp   = f2 + (3./XTemp).*f1 + (3./XTemp.^2).*f;
    Psi_tTemp = f3 + (3./XTemp).*f2 + (3./XTemp.^2).*f1;
    Psi_rTemp = -f3 - (3./XTemp.^2).*f1 - (3./XTemp).*f2 - (3./XTemp.^2).*f1 - (6./XTemp.^3).*f;
    
    temp = max(abs(PsiTemp(:))) + max(abs(Psi_tTemp(:))) + max(abs(Psi_rTemp(:)));
    %disp(['Maximum Violation of compact support = ',num2str(temp)])
    pause(2)
    
end