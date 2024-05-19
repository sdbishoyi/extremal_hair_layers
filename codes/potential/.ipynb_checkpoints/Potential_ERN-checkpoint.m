%rstar_ern = linspace(-100,1200,1300*5);
rstar_ern = -100;
M=0.5;
Q=M;
r_ern = RstarToR_ERN(rstar_ern,M);
Delta = (r_ern - M).^2 ;
pot = (2*M - 2*Q^2./r_ern).*(Delta./r_ern.^5);
pot(1)
%plot(r_rn, pot, 'LineWidth',1)