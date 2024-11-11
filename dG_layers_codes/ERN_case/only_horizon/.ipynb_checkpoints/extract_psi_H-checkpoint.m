psi_extr={};
extr_posn=s;
[r,c]=find(abs(x-extr_posn)<=1);
for i=1:length(tarr)
   psi_t=cell2mat(psiarr(i));
   psi=psi_t(r(1),c(1));
   psi_extr=[psi_extr,psi];

end
psi_extr=cell2mat(psi_extr);
% x_100=x;
h=figure;
plot(log10(tarr(4:end)),log10(abs(psi_extr(4:end))));
hold on
plot(log10(tarr(4:end)),-log10(tarr(4:end)))
hold off
%ylim([-15,1])
% plot properties
ylabel("log_{10} |\Psi (t,r*)|");
xlabel("log_{10}(t)")
t=sprintf("N=%d, subd=%d, FinalTime=%d, signal extracted at rstar=%1.4f",N,K,FinalTime,x(r(1),c(1)));
title(t)
%loc=strcat('plots/',t,'.png');
set(h,'Name' ,t);
%saveas(gcf,'Psivst_no_support_on_H_at_H.png')
saveas(gcf,'Psivst_support_on_Hinrho_at_H.png')
