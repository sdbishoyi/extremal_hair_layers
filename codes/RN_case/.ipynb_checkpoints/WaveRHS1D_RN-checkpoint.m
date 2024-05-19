function [rhsPsi,rhsPi, rhsPhi] = WaveRHS1D_RN(Psi,Pi,Phi,dg_globals,phys_system,time)



Np=dg_globals.Np;
Nfp=dg_globals.Nfp;
Nfaces=dg_globals.Nfaces;
K=dg_globals.K;
x=dg_globals.x;
vmapM=dg_globals.vmapM;
vmapP=dg_globals.vmapP;
nx=dg_globals.nx;
rx=dg_globals.rx;
Dr=dg_globals.Dr;
LIFT=dg_globals.LIFT;
Fscale=dg_globals.Fscale;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pot=phys_system.pot;
omega=phys_system.omega;
ohm=phys_system.ohm;
capH=phys_system.capH;
capHP=phys_system.capHP;
Potential=phys_system.Potential;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define field differences at faces
dPi = zeros(Nfp*Nfaces,K); dPi(:) = Pi(vmapM)-Pi(vmapP);
dPhi = zeros(Nfp*Nfaces,K); dPhi(:) = Phi(vmapM)-Phi(vmapP);


% %%Left boundary conditions
PiL= (Pi(1,1) -Phi(1,1))/2; dPi (1,1) = Pi(1,1) - PiL; 
PhiL =  (-Pi(1,1) +Phi(1,1))/2; dPhi (1,1) = Phi(1,1) - PhiL;


%%%Right boundary conditions
%PiR= (Pi(end,end) + Phi(end,end))/2; 
%dPi (end,end) = Pi(end,end) - PiR; 

%PhiR =  (Pi(end,end) +Phi(end,end))/2; 
%dPhi (end,end) = Phi(end,end) - PhiR;

%%%%%% evaluate upwind fluxes
% fluxPi = 0.5.*(nx.*dPhi - dPi);
% fluxPhi = 0.5.*(nx.*dPi - dPhi);



%%%%% calculate LF fluxes

Lambda=1;

if(phys_system.first_order_reduction_type==1)
    fluxPi = (1./(1+capH)).*(2.*capH.*Pi   + (1-capH).*Phi);
    fluxPhi= Pi;
elseif(phys_system.first_order_reduction_type==2)
    fluxPi = -1./(1+capH).*( - Phi - capH.*(Pi) ) ;
    fluxPhi =-1./(1+capH).*( -Pi - capH.*Phi);
else
    disp("ERROR! system type not found")
end

fluxPiM = fluxPi(vmapM);
fluxPiP = fluxPi(vmapP);


fluxPhiM = fluxPhi(vmapM);
fluxPhiP = fluxPhi(vmapP);

fluxPiP(end,end) = fluxPiM(end,end); 
fluxPhiP(end,end) = fluxPhiM(end,end); 

%% central numerical flux
fluxPitrue =(fluxPiP  + fluxPiM)/2;
fluxPhitrue=(fluxPhiP + fluxPhiM)/2;

%% Lax-Friedrich flux
fluxPitrue(:) = fluxPitrue(:)  +   (Lambda/2)*nx(:).*dPi(:);
fluxPhitrue(:) = fluxPhitrue(:)  +   (Lambda/2)*nx(:).*dPhi(:);

%%% NOTE TO MANAS: "Lax-Friedrich flux" - "F" 
%%% this is F*-F in the paper (See eq 35)
fluxPi = zeros(Nfp*Nfaces,K);   fluxPi(:) = fluxPitrue(:)   - fluxPiM(:) ;
fluxPhi = zeros(Nfp*Nfaces,K);  fluxPhi(:)= fluxPhitrue(:)  - fluxPhiM(:);

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


%%% NOTE TO MANAS: This original block of code commented out also works, 
%%% but the new block of code (below) is more compact and intuitive
%%% source term modifying subdomain to the left of Dirac delta
%%% Pi_source(2,rp_col_idx_left)=  (-G_00/2)*nx(2,rp_col_idx_left);
%%% Phi_source(2,rp_col_idx_left)= (G_00/2)*nx(2,rp_col_idx_left);

%%% source term modifying subdomain to the right of Dirac delta
%%% Pi_source(1,rp_col_idx_right)=  (G_00/2)*nx(1,rp_col_idx_right);
%%% Phi_source(1,rp_col_idx_right)= (G_00/2)*nx(1,rp_col_idx_right);

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
%%% modify the fluxes due to the Dirac delta
fluxPi =  -nx.*(fluxPi-Pi_source);
fluxPhi = -nx.*(fluxPhi-Phi_source);



% compute right hand sides of the PDEs
rhsPsi = -Pi;


% rhsPi = 1./(1+capH).*(   -(1-capH).*rx.*(Dr*Phi)...                             
%                             -2.*capH.*rx.*(Dr*Pi) ...
%                             + rx.*(Dr*capH).*Phi ...
%                             - rx.*(Dr*capH).*Pi) + LIFT*(Fscale.*fluxPi) ...  
%         - Potential.*Psi ;

if(phys_system.first_order_reduction_type==1)
    rhsPi = 1./(1+capH).*(   -(1-capH).*rx.*(Dr*Phi)...                             
                                -2.*capH.*rx.*(Dr*Pi) ...
                                + capHP.*Phi ...
                                - capHP.*Pi) + LIFT*(Fscale.*fluxPi) ...  
            - Potential.*Psi ;
    rhsPhi = (-rx.*(Dr*Pi) + LIFT*(Fscale.*fluxPhi));
elseif(phys_system.first_order_reduction_type==2)
    rhsPi = 1./(1+capH).*( - rx.*(Dr*Phi) - capH.*rx.*(Dr*Pi)) + LIFT*(Fscale.*fluxPi)...
        - Potential.*Psi ;
    rhsPhi =1./(1+capH).*( - rx.*(Dr*Pi) - capH.*rx.*(Dr*Phi)) + LIFT*(Fscale.*fluxPhi) ...
        -  capH.*Potential.*Psi ;
end


return

