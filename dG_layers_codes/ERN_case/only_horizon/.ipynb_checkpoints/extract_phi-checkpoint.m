phi_extr={};
extr_posn=800;
[r,c]=find(abs(x-extr_posn)<=1);
for i=1:length(tarr)
   phi_t=cell2mat(phiarr(i));
   phi=phi_t(r(1),c(1));
   phi_extr=[phi_extr,phi];
%    hold on
end
phi_extr=cell2mat(phi_extr);
% x_100=x;
h=figure;
plot(log10(tarr(4:end)),log10(abs(phi_extr(4:end))));
ylim([-50,5])
% plot properties
ylabel("log_{10} |\Phi (t,r*)|");
xlabel("log_{10}(t)")
t=sprintf("N=%d, subd=%d, FinalTime=%d, signal extracted at rstar=%1.4f",N,K,FinalTime,x(r(1),c(1)));
title(t)
%loc=strcat('plots/',t,'.png');
set(h,'Name' ,t);
%saveas(gcf,'Phivst_no_support_on_H_at_H.png')
saveas(gcf,'Phivst_support_on_Hinrho_at_rs=800.png')
