% The following is some code written by Manas for the evaluation of Potentials 
% Removed from : WaveDriver_ERN file. Added on : 1/14/24.

% different way of writting the potential, good for large r evaluations
% Pot_inf=((omega-x_64.*omegaP)./(1+capH))...
%         .*((x-2*M*omega)./x)...
%         .*((-2*M*omega)./x.^3 - ell.*(ell+1)./x.^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{ 
 Pot_inf = ((omega-x_64.*omegaP)./(1+capH))...
     .*((x.^2-2*M*omega.*x + M^2.*omega.^2)./x.^2)...
     .*(-(2*M - 2*M^2.*omega./x).*omega./x.^3 - ell.*(ell+1)./x.^2);
Potential_eff(Np,K)=double(Pot_inf(Np,K));
%}
% Potential_eff=-(omega-x.*omegaP).*(ell*(ell+1))./(x.*x.*(1+capH));
%Potential(end,end) = 0;
%Potential(end,end) = - Pot_large_r(end,end);

%%%%%%%%%%%% Some lines for the conversion routine. %%%%%%%%%%%%%
% Removed from : WaveDriver_ERN. Added on : 1/14/2024

%r_ern_som=RstarToRwithCharge_Newton(rstar_adj,M);
%r_ern(Np,K)=r_ern_som(Np,K);


%%%%%%%%%%%%%%%%%%% Some lines for the source term %%%%%%%%%%%%%%%%%%%%%%%%
% Removed from : WaveRHS1D_ERN . Added on 1/14/2024

%{
%%%%%%%%%% add particle
rp_rsch_value_final = phys_system.rp_rsch_value_final;
rp_row_idx=phys_system.rp_row_idx;
rp_col_idx_left=phys_system.rp_col_idx_left;
rp_col_idx_right=phys_system.rp_col_idx_right;

a=0;
M=1;
theta_p = pi/2;
vel=sqrt(M/rp_rsch_value_final);
atild=a/M;
omega_phi= vel^3/(M*(1+atild*vel^3));
phi_p = omega_phi*time;

mass_particle = 0.0;

[G00src] = SourceCoeff_Kerr_Teukolsky_s0(phys_system.ell, phys_system.m, theta_p, phi_p,rp_rsch_value_final, M, a, mass_particle);
smoother_end=400; smoother_delta=.00025;

%G00src=0.1*sin(0.1*time);

[mu, mu_t , mu_tt] = JumpSmoother(time,0,smoother_end,smoother_delta);
% disp(mu)
%mu=1;
G_00=1*mu*G00src;

% sourcevec=[G_00 0];

Pi_source = zeros(2,K);
Phi_source = zeros(2,K);

%%% NOTE TO MANAS: 
%%% fluxPi and fluxPi are components of F*-F (see above)
%%% Eq 37 (current version of the paper) shows how F* should be modified
%%% when we have non-zero source terms. 
%%%
%%% The coefficients (\pm 1/2) used below arise from acting on source vector [G_00,0] with various matricies
%%% Please see doc/Hyperbolicity-And-Source-TermJump.nb for how this is done


%%% modify the fluxes due to the Dirac delta
%%% fluxPi = -nx.*(fluxPi)+Pi_source;
%%% fluxPhi = -nx.*(fluxPhi)+Phi_source;


% source term modifying subdomain to the left of Dirac delta
Pi_source(2,rp_col_idx_left)=  -G_00/2;
Phi_source(2,rp_col_idx_left)=  G_00/2;

% source term modifying subdomain to the right of Dirac delta
Pi_source(1,rp_col_idx_right)=  G_00/2;
Phi_source(1,rp_col_idx_right)= G_00/2;

%%% NOTE TO MANAS: 
%%% The various negative signs can be understood as follows:
%%% 1. Note that fluxPi and fluxPhi are components of F*-F.
%%% 2. the "-" in front of Pi/Phi_source can be seen in Eq 37 (current draft)
%%%    and considering point 1 above
%%% 3. The negative sign in front of nx is because F-F* is the term that shows
%%% up on the right-hand-side of the PDF (see Eq 35 in current draft)
%%% 
%}

%%% modify the fluxes due to the Dirac delta
%fluxPi =  -nx.*(fluxPi-Pi_source);
%fluxPhi = -nx.*(fluxPhi-Phi_source);