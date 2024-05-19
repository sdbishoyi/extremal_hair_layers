init1 = 1500;

%%% We can try log10 or just log... lets make things systematic.
logphi1 = log10(abs(phi_extr(init1:end)));
%psi_extr_tail = psi_extr(init1:end);

%%% The following are the time arrays.
logt1 = log10(tarr(init1:end));
%tarr_tail = tarr(init1:end)-s;
%logt1 = log10(abs(tarr(init1:end)-s)/M);


%%% Now, we calculate the slope

%%% This is dlog(psi)dlog(tau).
slopephi = (logphi1(2:end)-logphi1(1:end-1))./(logt1(2:end)-logt1(1:end-1));
%%% This is tau.\dot{psi}/psi.
%{
slope2 = (tarr_tail(2:end)).*(psi_extr_tail(2:end) - psi_extr_tail(1:end-1))...
    ./(tarr_tail(2:end) - tarr_tail(1:end-1)).*1./(psi_extr_tail(2:end));

% logt1 = logt1(2:end);
% invt1 = 1./logt1;

%%% This is the ratio of the coefficients C2/C1.
coeff1 = (slope1+2).*(tarr(init1+1:end)-s);
coeff2 = -(slope1+2).*(tarr(init1+1:end)).*(log(tarr(init1+1:end))).^-1;
coeff3 = (-1/(4*M))*(coeff2+4*M);

%%% Now we plot.
%hold on
%plot((tarr(init1+1:end))/M,coeff3,"Color",'r','LineStyle','-','LineWidth',1)
plot(log(tarr(init1+1:end)-s),coeff1,"Color",'m','LineStyle','-')
%}

tail=figure;
plot(log10(tarr(init1+1:end)),slopephi,"Color",'b','LineStyle','-')
hold on
plot(log10(tarr(init1+1:end)),0*ones(length(tarr(init1+1:end))),"Color",'r','LineStyle','-')

%ylim([-5.0,-1.95])
%yticks([-2.10 -2.05 -2.00 -1.95 -1.90])
%ticks = linspace(-5.0,-1.95,10);
%yticks(ticks)
ylabel('LPI of \Phi')
%ylabel("C2/C1");
xlabel("log(u)")
tit=sprintf("N=%d, subd=%d, L=%d, FinalTime=%d, extracted at r*=%1.4f",N,K,ell,FinalTime,x(r(1),c(1)));
title(tit)
set(tail,'Name' ,tit);
%loc=strcat('plots/',t,'.png');
hold off
%saveas(gcf,'LPIvst_no_support_on_H_at_H.png')
saveas(gcf,'LPIphivst_support_on_Hinrho_at_rs=800.png')