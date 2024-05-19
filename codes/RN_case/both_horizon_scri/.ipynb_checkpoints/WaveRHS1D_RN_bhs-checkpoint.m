function [rhsPsi,rhsPi, rhsPhi] = WaveRHS1D_RN_bhs(Psi,Pi,Phi,dg_globals,phys_system,time)

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
F = dg_globals.F; % <-- notice the change here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pot=phys_system.pot;
omega=phys_system.omega;
ohm=phys_system.ohm;
capH=phys_system.capH;
% capHP=phys_system.capHP;
Potential=phys_system.Potential;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define field differences at faces
dPi = zeros(Nfp*Nfaces,K); dPi(:) = Pi(vmapM)-Pi(vmapP);
dPhi = zeros(Nfp*Nfaces,K); dPhi(:) = Phi(vmapM)-Phi(vmapP);


% %%Left boundary conditions
PiL= (Pi(1,1) -Phi(1,1))/2; dPi(1,1) = Pi(1,1) - PiL; 
PhiL =  (-Pi(1,1) +Phi(1,1))/2; dPhi (1,1) = Phi(1,1) - PhiL;


%%%Right boundary conditions
PiR= (Pi(end,end) + Phi(end,end))/2;  dPi (end,end) = Pi(end,end) - PiR; 
PhiR =  (Pi(end,end) +Phi(end,end))/2; dPhi (end,end) = Phi(end,end) - PhiR;


%%%%%% evaluate upwind fluxes
% fluxPi = 0.5.*(nx.*dPhi - dPi); fluxPhi = 0.5.*(nx.*dPi - dPhi);

%%%%% calculate LF fluxes

Lambda=1;

if(phys_system.first_order_reduction_type==1)
    fluxPi = (1./(1+capH)).*(2.*capH.*Pi   + (1-capH).*Phi);
    fluxPhi= Pi;
    
%changes for both layers done below on 2/27/2024. 
    
elseif(phys_system.first_order_reduction_type==2)
    fluxPi = -1./(1+capH).*( - Phi - capH.*F.*(Pi) ) ;
    fluxPhi =-1./(1+capH).*( -Pi - capH.*F.*Phi);
else
    disp("ERROR! system type not found")
end

fluxPiM = fluxPi(vmapM);
fluxPiP = fluxPi(vmapP);

fluxPhiM = fluxPhi(vmapM);
fluxPhiP = fluxPhi(vmapP);

fluxPiP(end,end) = fluxPiM(end,end); 
fluxPhiP(end,end) = fluxPhiM(end,end); 

fluxPiP(1,1) = fluxPiM(1,1);
fluxPhiP(1,1) = fluxPhiM(1,1);

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

fluxPi =  -nx.*(fluxPi);
fluxPhi = -nx.*(fluxPhi);

% compute right hand sides of the PDEs
rhsPsi = -Pi;

if(phys_system.first_order_reduction_type==1)
    rhsPi = 1./(1+capH).*(   -(1-capH).*rx.*(Dr*Phi)...                             
                                -2.*capH.*rx.*(Dr*Pi) ...
                                + capHP.*Phi ...
                                - capHP.*Pi) + LIFT*(Fscale.*fluxPi) ...  
            - Potential.*Psi ;
    rhsPhi = (-rx.*(Dr*Pi) + LIFT*(Fscale.*fluxPhi));
    
% changes for both layer done below on 2/27/2024.
    
elseif(phys_system.first_order_reduction_type==2)
    rhsPi = 1./(1+capH).*( - rx.*(Dr*Phi) - capH.*F.*rx.*(Dr*Pi)) + LIFT*(Fscale.*fluxPi)...
        - Potential.*Psi ;
    rhsPhi =1./(1+capH).*( - rx.*(Dr*Pi) - capH.*F.*rx.*(Dr*Phi)) + LIFT*(Fscale.*fluxPhi) ...
        -  capH.*F.*Potential.*Psi ;
end


return

