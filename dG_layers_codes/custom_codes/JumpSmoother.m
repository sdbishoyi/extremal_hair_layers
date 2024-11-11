function [mu, mu_t, mu_tt] = JumpSmoother(t,t0,tau,delta)

% used to smoothly turn on distributional source terms

if t < t0
  mu   =  0.0;
  mu_t =  0.0;
  mu_tt = 0.0;
else
  if t >= t0 + tau
     mu   =  1.0;
     mu_t =  0.0;
     mu_tt = 0.0;
  else
     mu   = 0.5*(erf(sqrt(delta)*(t-t0-0.5*tau))+1);
     mu_t = sqrt(delta/pi)*exp(-delta*(t-t0-0.5*tau).^2);
     mu_tt = -2*delta^(3/2)/sqrt(pi)*exp(-delta*(t-t0-.5*tau).^2).*(t-t0-.5*tau);
  end
end

%this funcion smoothly interpolates such that at t=0 f=f'=f''=f'''=f''''=0,
%and at t=1 f=1 f'=f''=f'''=f''''=0

% function [mu, mu_t, mu_tt] = JumpSmoother(t)
% if t <= 100
%     x=t/100;
%     mu = x.^5-5*x.^5.*(x-1)+15*x.^5.*(x-1).^2-35*x.^5.*(x-1).^3+70.*x.^5.*(x-1).^4;
%     mu_t = 5*x.^4+...
%     -25*x.^4.*(x-1)-5*x.^5 +...
%     5*15*x.^4.*(x-1).^2+2*15*x.^5.*(x-1)+...
%     -5*35*x.^4.*(x-1).^3-3*35*x.^5.*(x-1).^2+...
%     5*70.*x.^4.*(x-1).^4+4*70.*x.^5.*(x-1).^3;
%     mu_tt=0;
% else
%      mu   =  1.0;
%      mu_t =  0.0;
%      mu_tt = 0.0;
% end
