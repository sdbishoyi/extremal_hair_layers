% This file contains the potential for sub-extremal RN BHs.
% The formula can be noted from : eqn(4). https://arxiv.org/pdf/1707.00515.pdf

%rstar_rn = linspace(-100,1200,1300*5);
rstar_rn = -100;
M=0.5;
Q=0.3;
r_rn = RstarToR_RN(rstar_rn,M,Q);
Delta = r_rn.^2 -2*M.*r_rn + Q^2;
pot_rn = (l(l+1)./r_rn + 2*M./r_rn - 2*Q^2./r_rn).*(Delta./r_rn.^4);
pot_rn(1)
%plot(r_rn, pot_rn, 'LineWidth',1)