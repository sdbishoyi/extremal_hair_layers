init1 = 5000;

%%% We can try log10 or just log... lets make things systematic.
logpsi1 = log(abs(psi_extr(init1:end)));
psi_extr_tail = psi_extr(init1:end);

%%% The following are the time arrays.
logt1 = log(tarr(init1:end)-s);
tarr_tail = tarr(init1:end)-s;
%logt1 = log10(abs(tarr(init1:end)-s)/M);


%%% Now, we calculate the slope
%%% This is dlog(psi)dlog(tau).
slope1 = (logpsi1(2:end)-logpsi1(1:end-1))./(logt1(2:end)-logt1(1:end-1));
%%% This is tau.\dot{psi}/psi.
slope2 = (tarr_tail(2:end)).*(psi_extr_tail(2:end) - psi_extr_tail(1:end-1))...
    ./(tarr_tail(2:end) - tarr_tail(1:end-1)).*1./(psi_extr_tail(2:end));

% logt1 = logt1(2:end);
% invt1 = 1./logt1;

%%% This is the ratio of the coefficients C2/C1.
coeff = (slope1+2).*(tarr(init1+1:end)-s);
%coeffratio = -(slope1+2).*(tarr(init1+1:end)-s)./(log(tarr(init1+1:end)-s));
%%% Now we plot.
%hold on

%plot(1./log(tarr(init1+1:end)-s),coeffratio,"Color",'r','LineStyle','-','LineWidth',1)
plot(log(tarr(init1+1:end)-s),coeff,"Color",'m','LineStyle','-')
%plot(log(tarr(init1+1:end)-s),slope1,"Color",'m','LineStyle','-')
ylabel("C2/C1");
xlabel("log(u)")
% tit=sprintf("N=%d, subd=%d, L=%d, FinalTime=%d, extracted at r*=%1.4f",N,subd,ell,FinalTime,x(r(1),c(1)));
% title(tit)
% hold off