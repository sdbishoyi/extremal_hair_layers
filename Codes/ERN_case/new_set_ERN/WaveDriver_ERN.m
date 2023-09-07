%addpath(genpath("../../../00_essential_codes/"));

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
M=1;
Q=M;
a=0;
N=12;

FinalTime = 10000;

ell = 0;
m=0;
%xL = -100;
%xR = 200;
disp(["The mode being solved is ",num2str(ell),",",num2str(m)])

%%%% Make DG grid around the particle
rp_r=10.0;
rh_plus=M + sqrt(M^2 - a^2);
rh_minus=M - sqrt(M^2 - a^2); 
rp_rstar=  rp_r + (2*M*rh_plus/(rh_plus-rh_minus))*log(abs(rp_r-rh_plus)/2*M) ...
            - (2*M*rh_minus/(rh_plus-rh_minus))*log(abs(rp_r-rh_minus)/2*M);
            
%rp_rstar should be an interface in the grid.

xL=rp_rstar-200;
xR=rp_rstar+800;
x_regions=[xL rp_rstar xR];  %%% make chunks of dG subdomains

k1=100;
k2=400;

K_regions=[k1 k2]; % subdomains to use to the left/right of rp
[Nv, VX, K, EToV] = MeshGen1D_V1(x_regions,K_regions);
StartUp1D;

xmin_left_rp= min(diff(x(:,1)));
xmin_right_rp= min(diff(x(:,1+k1)));

assert(abs(xmin_right_rp-xmin_left_rp)<=1e-10)
dg_globals.xmin=min([xmin_left_rp xmin_right_rp]);

%%%%


%%%% Generate simple mesh

%subd=200; 

%[Nv, VX, K, EToV] = MeshGen1D(xL,xR,subd);
%%%%
% Initialize solver and construct grid and metric
% StartUp1D;
% fid = fopen('./x_values_matlab.dat','w');
% fprintf(fid,'%32.32f \n',x); 
% fclose(fid);

%locR=30;
locR=200;


[idx1,idx2]=find_value_arr(x,locR);
%[idx1,idx2]=find(abs(x-locR)<1);
idx1=idx1(1); idx2=idx2(1);
locR=x(idx1,idx2);
s=xR;
P=4;


omega = 1  - ( ( (x-locR)./(s-locR) ).^P).*homeHVSD(x-locR);

omegaP = (-1).*(((-1).*locR+s).^(-1).*((-1).*locR+x)).^P.*0 ...
                +(-1).*P.*((-1).*locR+s)^(-1).*(((-1).*locR+s).^(-1) ...
                *((-1).*locR+x)).^((-1)+P).*homeHVSD((-1).*locR+x);
ohm = omega.^2./(omega - x.*omegaP);

capH   = 1-ohm;

% capHP=-ohmP;

dg_globals.x=x;           dg_globals.Np=Np;        dg_globals.Nfp=Nfp; 
dg_globals.Nfaces=Nfaces; dg_globals.K=K;          dg_globals.vmapM=vmapM; 
dg_globals.vmapP=vmapP;   dg_globals.nx=nx;        dg_globals.rx=rx;
dg_globals.Dr=Dr;         dg_globals.LIFT=LIFT;    dg_globals.Fscale=Fscale;
dg_globals.rk4a=rk4a;     dg_globals.rk4b=rk4b;    dg_globals.rk4c=rk4c;



%%%%%%% 64 digits calculation
digits(32)
locR_64 = vpa(locR);
s_64 = vpa(s); 
P_64 = vpa(P);
x_64 = vpa(x);
C=1;
omega_64 = 1  - ( ( (x_64-locR_64)./(s_64-locR_64) ).^P_64).*homeHVSD(x_64-locR_64);


omegaP_64 = (-1).*(((-1).*locR_64+s_64).^(-1).*((-1).*locR_64+x_64)).^P_64.*0 ...
                +(-1).*P_64.*((-1).*locR_64+s_64)^(-1).*(((-1).*locR_64+s_64).^(-1) ...
                *((-1).*locR_64+x_64)).^((-1)+P_64).*homeHVSD((-1).*locR_64+x_64);
ohm_64 = omega_64.^2./(omega_64 - x_64.*omegaP_64);

ohmP_64 =x_64.*((-2).*C.*P_64.*((-1).*locR_64+s_64).^(-1).*(((-1).*locR_64+s_64).^(-1).*((-1) ...
  .*locR_64+x_64)).^((-1)+P_64).*0+(-1).*C.*((-1)+P_64).* ...
  P_64.*((-1).*locR_64+s_64).^(-2).*(((-1).*locR_64+s_64).^(-1).*((-1).*locR_64+x_64)).^( ...
  (-2)+P_64).*homeHVSD((-1).*locR_64+x_64)).*(1+(-1).*C.*(((-1).*locR_64+s_64).^( ...
  -1).*((-1).*locR_64+x_64)).^P_64.*homeHVSD((-1).*locR_64+x_64)).^2.*(1+(-1).*C.*( ...
  ((-1).*locR_64+s_64).^(-1).*((-1).*locR_64+x_64)).^P_64.*homeHVSD((-1).*locR_64+x_64)+( ...
  -1).*x_64.*((-1).*C.*(((-1).*locR_64+s_64).^(-1).*((-1).*locR_64+x_64)).^P_64.* ...
  0+(-1).*C.*P_64.*((-1).*locR_64+s_64).^(-1).*(((-1).* ...
  locR_64+s_64).^(-1).*((-1).*locR_64+x_64)).^((-1)+P_64).*homeHVSD((-1).*locR_64+x_64))) ...
  .^(-2)+2.*((-1).*C.*(((-1).*locR_64+s_64).^(-1).*((-1).*locR_64+x_64)).^P_64.* ...
  0+(-1).*C.*P_64.*((-1).*locR_64+s_64).^(-1).*(((-1).* ...
  locR_64+s_64).^(-1).*((-1).*locR_64+x_64)).^((-1)+P_64).*homeHVSD((-1).*locR_64+x_64)) ...
  .*(1+(-1).*C.*(((-1).*locR_64+s_64).^(-1).*((-1).*locR_64+x_64)).^P_64.*homeHVSD( ...
  (-1).*locR_64+x_64)).*(1+(-1).*C.*(((-1).*locR_64+s_64).^(-1).*((-1).*locR_64+x_64)) ...
  .^P_64.*homeHVSD((-1).*locR_64+x_64)+(-1).*x_64.*((-1).*C.*(((-1).*locR_64+s_64).^( ...
  -1).*((-1).*locR_64+x_64)).^P_64.*0+(-1).*C.*P_64.*((-1) ...
  .*locR_64+s_64).^(-1).*(((-1).*locR_64+s_64).^(-1).*((-1).*locR_64+x_64)).^((-1)+P_64) ...
  .*homeHVSD((-1).*locR_64+x_64))).^(-1);

capH_64   = 1-ohm_64;

capHP_64=-ohmP_64;


r_adj_64=x_64./omega_64;

harr_64=(r_adj_64./omega_64) - r_adj_64;
% pot_64 = -ell*(ell+1)./r_adj_64.^2;

%%%%

omega=double(omega_64);
ohm=double(ohm_64);
ohmP=double(ohmP_64);
omegaP=double(omegaP_64);
capH=double(capH_64);
capHP=double(capHP_64);
r_adj=double(r_adj_64);
harr=double(harr_64);
% pot=double(pot_64);

rstar_adj = x./(omega); % need to prevent Inf at x(end,end)
harr=(rstar_adj./omega) - rstar_adj;

% pot = -ell*(ell+1)./r_adj.^2;

%[r_sch,r_schm2m] = RstarToR(rstar_adj,M);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[r_sch,r_schm2m] = RstarToR_ERN(rstar_adj,M);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r_sch_manas=RstarToRwithSpin_Newton(rstar_adj,M,a);
r_sch(Np,K)=r_sch_manas(Np,K);
%pot_sch=(1-2*M./r_sch).*(-2*M./(r_sch.^3) - ell*(ell+1)./r_sch.^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pot_sch=-(ell*(ell + 1)./r_sch.^4 + (2*M - 2*M^2./r_sch)./r_sch.^5).*...
     (r_sch.^2 -2*M.*r_sch + M^2);

% mathematicavalues;
% omega=omega_mathematica;
% omegaP=omegaP_mathematica;
% ohm=ohm_mathematica;
% r_adj=rstar_mathematica;
% ohmP=ohmP_mathematica;
% capH=capH_mathematica;
% capHP=capHP_mathematica;


Potential_eff = (1./((1-capH).*(1+capH))).*pot_sch;


% Potential=Potential_mathematica;
% different way of writting the potential, good for large r evaluations

%%% Manas commented this out as it uses r_adj
% Pot_large_r = (1./r_adj - omegaP).*(ell*(ell+1))./(x.*(1+capH)

% Pot_inf=((omega-x_64.*omegaP)./(1+capH))...
%         .*((x-2*M*omega)./x)...
%         .*((-2*M*omega)./x.^3 - ell.*(ell+1)./x.^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 Pot_inf = ((omega-x_64.*omegaP)./(1+capH))...
     .*((x.^2-2*M*omega.*x + M^2.*omega.^2)./x.^2)...
     .*(-(2*M - 2*M^2.*omega./x).*omega./x.^3 - ell.*(ell+1)./x.^2);


Potential_eff(Np,K)=double(Pot_inf(Np,K));

% Potential_eff=-(omega-x.*omegaP).*(ell*(ell+1))./(x.*x.*(1+capH));

%Potential(end,end) = 0;
%Potential(end,end) = - Pot_large_r(end,end);
Potential=Potential_eff;
phys_system.Potential=Potential;


%%%%%%%%%%%%%% Set initial condition for actual problem %%%%%%%%%%%%%%%%%%%

sigma=9; muu=3;
%%% Below is static IC
% psi_in = (1/sqrt(2*pi*sigma^2))*exp(-(x-muu).^2/(2*sigma^2));
% pi_in = zeros(Np,K);
% phi_in = (1/sqrt(2*pi*sigma^2))*exp(-(x-muu).^2/(2*sigma^2)).*(-2*(x-muu)./(2*sigma^2)); 

%%% Below is generic IC type 1
% psi_in = 0*(1/sqrt(2*pi*sigma^2))*exp(-(x-muu).^2/(2*sigma^2));
% pi_in = (1/sqrt(2*pi*sigma^2))*exp(-(x-muu).^2/(2*sigma^2));
% phi_in = 0*(1/sqrt(2*pi*sigma^2))*exp(-(x-muu).^2/(2*sigma^2)).*(-2*(x-muu)./(2*sigma^2)); 

%%% Below is generic IC type 2
psi_in = 1*(1/sqrt(2*pi*sigma^2))*exp(-(x-muu).^2/(2*sigma^2));
pi_in = -1*(1/sqrt(2*pi*sigma^2))*exp(-(x-muu).^2/(2*sigma^2)).*(-2*(x-muu)./(2*sigma^2));
phi_in = -pi_in;

%%%%%% SCOTT: exact outgoing initial data for ell=2
%x0=-10;
%f0=2;
%c=1;
%[Psi,Psi_r,Psi_t,Psi_tt] = OutgoingEll2_ID_SineGaussian(x,f0,c,x0,-x);
%psi_in = Psi;
%pi_in = -Psi_t;
%phi_in = Psi_r;


%%%%%%%zero IC
% psi_in=zeros(size(x));
% pi_in=zeros(size(x));
% phi_in=zeros(size(x));


init_con.psi=psi_in;     init_con.pi=pi_in;        init_con.phi=phi_in;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   physical system%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phys_system.ell=ell;        
% phys_system.pot=pot;
phys_system.omega=omega;
phys_system.omegaP=omegaP;
phys_system.ohm=ohm;
phys_system.capH=capH;
phys_system.capHP=capHP;
phys_system.first_order_reduction_type = first_order_reduction_type;
phys_system.ell=ell;
phys_system.m=m;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set spectral filter
%Filter_Nc = floor(2*N/3); %cut on modes above this value
%Filter_s = 8; %exp filter, must be even number
%matF=Filter1D(dg_globals,Filter_Nc,Filter_s);
matF=1;
dg_globals.matF=matF;
dg_globals.r_adj=r_adj;






%%%%%%%%%%% particle at rp
% rp_rsch_value_init=10;
% [rp_row_idx,rp_col_idx]=find_value_arr(r_sch,rp_rsch_value_init);
% rp_row_idx=Np;
% rp_col_idx_left=rp_col_idx; rp_col_idx_right=rp_col_idx+1; 
rp_row_idx=Np; 
rp_col_idx_left=k1;
rp_col_idx_right=k1+1;


rp_rsch_value_final=r_sch(rp_row_idx,rp_col_idx_left);
disp(["Particle is at r = ", num2str(rp_rsch_value_final,32) , " and r* = ",num2str(x(rp_row_idx,rp_col_idx_left),32)])

phys_system.rp_rsch_value_final=rp_rsch_value_final;
phys_system.rp_row_idx=rp_row_idx;
phys_system.rp_col_idx_left=rp_col_idx_left;
phys_system.rp_col_idx_right=rp_col_idx_right;


[psiarr,piarr,phiarr,tarr] = Wave1D_ERN(init_con,dg_globals,phys_system,FinalTime);
tarr=cell2mat(tarr);
psiarr=psiarr(~cellfun(@isempty, psiarr));

piarr=piarr(~cellfun(@isempty, piarr));
phiarr=phiarr(~cellfun(@isempty, phiarr));


