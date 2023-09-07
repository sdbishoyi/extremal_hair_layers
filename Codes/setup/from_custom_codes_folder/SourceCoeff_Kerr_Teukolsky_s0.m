function [G] = SourceCoeff_Kerr_Teukolsky_s0(ell, m, theta_p, phi_p, rp, mass_kerr, spin_kerr, mass_particle)

%
% for circular orbits, E, L found from spin and rp. Generic function below
% might be useful later on
% function [G] = SourceCoeff_Kerr_Teukolsky_s0(ell, m, rp, E_orbital, L_orbital, mass_kerr, spin_kerr, mass_particle)


% This function computes the G_{ell,m}(t) source coefficient that arises
% on the left-hand-side of the s=0 Teukolsky equation when the perturbing
% particle is geodesic in circular orbit located at \theta = pi/2
%
% Here the Teukolsky has been decomposed into coupled (ell,m) modes
% in Boyer-Lindquist coordinates.
%
% Scalar charge density T is eq 8 of https://arxiv.org/pdf/1003.1860.pdf
%
%  Input:
%     ell, m are the spherical indicies 
%     rp is the location of the particle in r
%     theta_p, phi_p is the location of the theta and phi
%     E_orbital, L_orbital are the orbit's energy and angular momentum 
%     mass_kerr, spin_kerr are parameters of the background Kerr blackhole 
%     mass_particle   is the mass of the orbiting test particle 
%
%
%

assert(spin_kerr<1.0)
assert(mass_kerr==1)  % equations assume 1
assert(mass_particle<=1)
assert(theta_p==pi/2)



%% Eq 2 of https://arxiv.org/pdf/1003.1860.pdf
v = sqrt(mass_kerr/rp);
a = spin_kerr;
E_orbital = (1-2*v^2+a*v^3)/sqrt(1-3*v^2+2*a*v^3);
L_orbital = rp*v*(1-2*a*v^3+a*a*v^4)/sqrt(1-3*v^2+2*a*v^3);

%% t-component of the 4-velocity vector
M=mass_kerr;
Deltarp=rp*rp + a*a - 2*M*rp;

%%%%%%% Kerr parameters
Sig=rp*rp + a*a*cos(theta_p)*cos(theta_p);
g_t_phi_up = -2*M*rp*a/(Sig*Deltarp*sin(theta_p)*sin(theta_p)); % TODO FOR MANAS... ADD CORRECT FACTOR % done, used  https://www.roma1.infn.it/teongrav/onde19_20/kerr.pdf as reference
g_t_t_up = -(1/Deltarp)*(rp*rp + a*a + (1/Sig)*2*M*rp*a*a*sin(theta_p)*sin(theta_p)); % TODO FOR MANAS... ADD CORRECT FACTOR


%%%%%% define ut at the location of particle
ut = g_t_phi_up*L_orbital - g_t_t_up*E_orbital;





ylm = SphericalHarmonics_HighEll(ell,m,theta_p,phi_p);
% G = - ( 1/(ut*sqrt(rp*rp + a*a)) )*(4*pi*mass_particle*(rp*rp + cos(theta_p)*cos(theta_p)))*conj(ylm);
G = - ( 1/(ut*sqrt(rp*rp + a*a)) )*(4*pi*mass_particle)*conj(ylm);

%%%% transform delta function to tortoisee coords 
% drs_dr=(rp*rp + a*a)/Deltarp;
% G=G*drs_dr;




% G=2*G/Deltarp;




% G=1*G;  %%% coefficient should be twice the jump
% G=G/(Deltarp);  %%%% Barry's suggestion because Niels' paper uses r*R for jump calculations
% drs_dr=(rp*rp + a*a)/Deltarp;
% G=G*drs_dr;
% G=-G;
% G=G/(4*pi);   %%%%   still not sure about this change...need to double check
end