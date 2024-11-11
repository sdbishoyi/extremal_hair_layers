addpath('../setup')
addpath('../custom_codes')
addpath('../conversion')
Globals1D;

% NOTE: WaveRHS1D impliments two different versions of the wave equation 
% using hlayers.
%
% System 1: uses phi = \partial_{rho} psi
%           uses pi  = \partial_{tau} psi
%
% System 2: uses phi = \partial_{r} psi
%           uses pi  = \partial_{t} psi
%
first_order_reduction_type=2;

% Polynomial order used for approximation 
FinalTime = 2000;
M=1; Q=M; a=0; ell = 0; 


N=9; xL = -10; xR = 1000; K = 600;

%disp(["The mode being solved is ",num2str(ell),",",num2str(m)])

[Nv, VX, K, EToV] = MeshGen1D(xL,xR,K);
StartUp1D;

xmin = min(abs(x(1,:)-x(2,:)))
dg_globals.xmin=xmin;
locR=-5;


[idx1,idx2]=find_value_arr(x,locR);
%[idx1,idx2]=find(abs(x-locR)<1);
idx1=idx1(1); idx2=idx2(1);
locR=x(idx1,idx2);
s=xL;
P=4;


omega = 1  - ( ( (x-locR)./(s-locR) ).^P).*homeHVSD(locR-x);

omegaP = (-1).*(((-1).*locR+s).^(-1).*((-1).*locR+x)).^P.*0 ...
                +(-1).*P.*((-1).*locR+s)^(-1).*(((-1).*locR+s).^(-1) ...
                *((-1).*locR+x)).^((-1)+P).*homeHVSD((+1).*locR-x);
                
ohm = omega.^2./(omega - x.*omegaP);

capH   = 1-ohm;

dg_globals.x=x;           dg_globals.Np=Np;        dg_globals.Nfp=Nfp; 
dg_globals.Nfaces=Nfaces; dg_globals.K=K;          dg_globals.vmapM=vmapM; 
dg_globals.vmapP=vmapP;   dg_globals.nx=nx;        dg_globals.rx=rx;
dg_globals.Dr=Dr;         dg_globals.LIFT=LIFT;    dg_globals.Fscale=Fscale;
dg_globals.rk4a=rk4a;     dg_globals.rk4b=rk4b;    dg_globals.rk4c=rk4c;

rstar_adj = x./(omega); % need to prevent Inf at x(end,end)

%%%%%%%%%%%%%%%% Set conversion routines and potentials %%%%%%%%%%%%%%%%

r_ern=RstarToRwithCharge_Newton(rstar_adj,M); % -- changed on 12/22/23

%pot_sch=(1-2*M./r_sch).*(-2*M./(r_sch.^3) - ell*(ell+1)./r_sch.^2);


pot_ern=-(ell*(ell + 1)./r_ern.^4 + (2*M - 2*M^2./r_ern)./r_ern.^5).*...
     (r_ern.^2 -2*M.*r_ern + M^2); %-- checked on 12/22/23

Potential_eff = (1./((1-capH).*(1+capH))).*pot_ern;
Potential_eff(1,1)=double(0);


Potential=Potential_eff;
phys_system.Potential=Potential;

%%%%%%%%%%%%%% Set initial condition for actual problem %%%%%%%%%%%%%%%%%%%


sigma=3.2605; muu=-2.7233 ;% --> in 'x' that is rho coordinate -- no support on H.

%%% Below is static IC
psi_in = (1/sqrt(2*pi*sigma^2))*exp(-(x-muu).^2/(2*sigma^2));
pi_in = zeros(Np,K);
phi_in = (1/sqrt(2*pi*sigma^2))*exp(-(x-muu).^2/(2*sigma^2)).*(-2*(x-muu)./(2*sigma^2)).*ohm; % --> the factor of 'ohm' is necessary.

%%% Below is generic IC type 1
%psi_in = (1/sqrt(2*pi*sigma^2))*exp(-(x-muu).^2/(2*sigma^2));
%pi_in = -(1/sqrt(2*pi*sigma^2))*exp(-(x-muu).^2/(2*sigma^2)).*(-2*(x-muu)./(2*sigma^2)); 
%phi_in = (1/sqrt(2*pi*sigma^2))*exp(-(x-muu).^2/(2*sigma^2)).*(-2*(x-muu)./(2*sigma^2)); % --> correction made on 1/2/2024 (Am).

%%% Below is generic IC type 2
%psi_in = 0*(1/sqrt(2*pi*sigma^2))*exp(-(x-muu).^2/(2*sigma^2));
%pi_in = (1/sqrt(2*pi*sigma^2))*exp(-(x-muu).^2/(2*sigma^2));
%phi_in = 0*-pi_in;

%%%%%%%%%%%%%%%%% Initial conditions in r_ern coordinates %%%%%%%%%%%%%%%%%

%sigma=0.16; muu=1.05;  %--> in r_ern coordinates for support on Horizon.
% psi_in = (1/sqrt(2*pi*sigma^2))*exp(-(r_ern-muu).^2/(2*sigma^2));
% pi_in = (1/sqrt(2*pi*sigma^2))*exp(-(r_ern-muu).^2/(2*sigma^2));
% pi_in = zeros(Np,K);
% phi_in = (1/sqrt(2*pi*sigma^2))*exp(-(r_ern-muu).^2/(2*sigma^2)).*(-2*(x-muu)./(2*sigma^2)); 

init_con.psi=psi_in;     init_con.pi=pi_in;        init_con.phi=phi_in;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% physical system %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phys_system.ell=ell;        
% phys_system.pot=pot;
phys_system.omega=omega;
phys_system.omegaP=omegaP;
phys_system.ohm=ohm;
phys_system.capH=capH;
%phys_system.capHP=capHP;
phys_system.first_order_reduction_type = first_order_reduction_type;
phys_system.ell=ell;
%phys_system.m=m;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set spectral filter
%Filter_Nc = floor(2*N/3); %cut on modes above this value
%Filter_s = 8; %exp filter, must be even number
%matF=Filter1D(dg_globals,Filter_Nc,Filter_s);
matF=1;
dg_globals.matF=matF;
dg_globals.rstar_adj=rstar_adj;

% Solve Problem

[psiarr,piarr,phiarr,tarr] = Wave1D_ERN(init_con,dg_globals,phys_system,FinalTime);

tarr=cell2mat(tarr);
psiarr=psiarr(~cellfun(@isempty, psiarr));
piarr=piarr(~cellfun(@isempty, piarr));
phiarr=phiarr(~cellfun(@isempty, phiarr));


