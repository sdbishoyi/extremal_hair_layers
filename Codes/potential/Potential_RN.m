%rstar_rn = linspace(-100,1200,1300*5);
rstar_rn = -100;
M=0.5;
Q=0.3;
r_rn = RstarToR_RN(rstar_rn,M,Q);
Delta = r_rn.^2 -2*M.*r_rn + Q^2 ;
pot = (2*M - 2*Q^2./r_rn).*(Delta./r_rn.^5);
pot(1)
%plot(r_rn, pot, 'LineWidth',1)