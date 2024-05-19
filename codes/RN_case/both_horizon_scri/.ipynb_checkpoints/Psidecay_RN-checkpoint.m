%function Psidecay_RN(begin_tail)
init1 = 1800; %begin_tail;

%%% We can try log10 or just log... lets make things systematic.

logpsi1 = log10(abs(psi_extr(init1:end)));

%psi_extr_tail = psi_extr(init1:end);
%%% The following are the time arrays.

logt1 = log10(tarr(init1:end));

%tarr_tail = tarr(init1:end)-s;
%logt1 = log10(abs(tarr(init1:end)-s)/M);

%%% Now, we calculate the slope
%%% This is dlog(psi)dlog(tau).

slope1 = (logpsi1(2:end)-logpsi1(1:end-1))./(logt1(2:end)-logt1(1:end-1));

%%% This is tau.\dot{psi}/psi.

tail=figure;
plot(log10(tarr(init1+1:end)),slope1,'-b','LineWidth',1.75)
hold on
plot(log10(tarr(init1+1:end)),-3*ones(length(tarr(init1+1:end))),'-r','LineWidth',1.75)

%ylim([-3,0])
%yticks([-3.10 -3.05 -3.00 -2.95 -2.90])

ylabel('LPI')

%ylabel("C2/C1");

xlabel("log_{10}(\tau)")
tit=sprintf("N=%d, subd=%d, L=%d, FinalTime=%d, extracted at rho=%1.4f",N,K,ell,FinalTime,x(r(end),c(end)));
title(tit)

%saveas(gcf,'LPIvst_no_support_on_H_rs=500.png')

filename = fullfile('../plots', 'LPIvst_no_support_on_Hinrho_rho=500.png');
%filename = fullfile('../plots', 'LPIvst_no_support_on_Hinrho_rho=xR.png');
saveas(gcf,filename);
set(tail,'Name' ,tit);
hold off
%end
