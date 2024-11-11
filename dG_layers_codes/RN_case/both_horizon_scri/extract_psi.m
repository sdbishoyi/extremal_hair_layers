%function extract_psi(extract_signal_at)
psi_extr={};
extr_posn=500;
%extr_posn=xR;    %extract_signal_at;
[r,c]=find(abs(x-extr_posn)<=1);
for i=1:length(tarr)
   psi_t=cell2mat(psiarr(i));
   psi=psi_t(r(end),c(end));
   psi_extr=[psi_extr,psi];
end
psi_extr=cell2mat(psi_extr);
h=figure;
plot(log10(tarr(4:end)),log10(abs(psi_extr(4:end))),'LineWidth',1.75);
ylim([-20,0])
% plot properties
ylabel("log_{10} |\Psi (\tau,\rho)|");
xlabel("log_{10}(\tau)")
t=sprintf("N=%d, subd=%d, FinalTime=%d, signal extracted at rho=%1.4f",N,K,FinalTime,x(r(end),c(end)));
title(t)
set(h,'Name' ,t);
%saveas(gcf,'Psivst_no_support_on_H_rs=500.png')
filename = fullfile('../plots', 'Psivst_no_support_on_H_rho=500.png');
%filename = fullfile('../plots', 'Psivst_no_support_on_H_rho=xR.png');
saveas(gcf,filename);
%end
