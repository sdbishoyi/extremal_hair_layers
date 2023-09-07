function [Psi_arr,Pi_arr,Phi_arr,t_arr] = Wave1D_ERN(init_con,dg_globals,phys_system,FinalTime) 

x=dg_globals.x;
Np=dg_globals.Np;
K=dg_globals.K;
rk4a=dg_globals.rk4a;
rk4b=dg_globals.rk4b;
rk4c=dg_globals.rk4c;
matF=dg_globals.matF;

time = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Runge-Kutta residual storage  
resPsi = zeros(Np,K); resPi = zeros(Np,K); resPhi = zeros(Np,K); 

% compute time step size
%xmin = min(abs(x(1,:)-x(2,:)));
xmin=dg_globals.xmin;
CFL=0.65;  

dt=0.01046059193835873271649639093539;
%dt=0.0064349;
%dt = CFL*xmin/5;
%dt=dt/(0.5*phys_system.m);
%dt=1e-4;
%CFL=dt/xmin
DT=0.5;
next_snapshot = DT;

disp(['xmin = ', num2str(xmin,8), ' dt = ', num2str(dt,8)])
% dt=1e-4
Nsteps = ceil(FinalTime/dt); dt = FinalTime/Nsteps;


Psi=init_con.psi; Pi=init_con.pi;  Phi=init_con.phi;

Psi_arr=cell(1,Nsteps);     Psi_arr{1}=Psi;
Pi_arr=cell(1,Nsteps);      Pi_arr{1}=Pi;
Phi_arr=cell(1,Nsteps);     Phi_arr{1}=Phi;
t_arr=cell(1,Nsteps);       t_arr{1}=time;

tic % used to estimte remaining simulation time

for tstep=1:Nsteps
   for INTRK = 1:5
      timelocal = time + rk4c(INTRK)*dt;
      [rhsPsi,rhsPi, rhsPhi] = WaveRHS1D_ERN(Psi,Pi,Phi,dg_globals,phys_system,timelocal);
       
      resPsi = rk4a(INTRK)*resPsi + dt*(rhsPsi);     
      resPi  = rk4a(INTRK)*resPi  + dt*(rhsPi);
      resPhi = rk4a(INTRK)*resPhi + dt*(rhsPhi);

      Psi = Psi + rk4b(INTRK)*resPsi;     
      Pi  = Pi  + rk4b(INTRK)*resPi;
      Phi = Phi + rk4b(INTRK)*resPhi;
      
      Psi=matF*Psi;
      Pi=matF*Pi;
      Phi=matF*Phi;
   
   end 
   % Increment time
   time = tstep*dt;
   abs(next_snapshot - time);

   if mod(time,DT)<dt/2 && tstep~=1
%    if abs(next_snapshot - time)<dt/2 && tstep~=1
       Psi_arr{tstep}=Psi;
       Pi_arr{tstep}=Pi;
       Phi_arr{tstep}=Phi;
  
       t_arr{tstep}=time;
       next_snapshot = next_snapshot + DT;
   end
   
   percmsg=20;
   printTee=floor(FinalTime/percmsg);
%    printTee=0.75;
   if mod(timelocal,printTee) <= dt && timelocal > 2*dt
        elapsed = toc/60;
        total_est = FinalTime*elapsed/timelocal;
        remain_est = total_est - elapsed;
        disp([num2str(timelocal*100/FinalTime),'% completed ',' | Program time: ',num2str(round(timelocal)) ,' | Elapsed clock(m) : ',num2str(elapsed),' | Estimated Time Remaining : ', num2str(remain_est) ])
        %disp(['-----------------------------------------'])
        %disp([' '])
        %disp(['Program time       : ',num2str(round(timelocal))])
        %disp(['Elapsed clock  (m) : ',num2str(elapsed)])
        %disp(['Time remaining (m) : ',num2str(remain_est)])
        %disp([' '])
        %disp(['-----------------------------------------'])

        %plot(x,Psi)
        %xlim([-10 30])
        %drawnow

   end
   

end
disp('Process completed succesfully')
return
