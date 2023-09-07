%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
psi_extr={};
extr_posn = s;
[r,c]=find(abs(x-extr_posn)<=1);
r=r(end); c=c(end);
 for i=1:length(tarr)
    psi_t=cell2mat(psiarr(i));
    psi=psi_t(r,c);
    psi_extr=[psi_extr,psi];
 end
 psi_extr= cell2mat(psi_extr);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
plot(tarr(4:end),log10(abs(psi_extr(4:end))),"LineWidth",1.25,...
     "Color",'r',LineStyle='-');
 ylim([-15,5])
 ylabel("log|\Psi (t,r*)|");
 xlabel("t")
 t=sprintf("N=%d, subd=%d, FinalTime=%d, signal extracted at ..." + ...
     "rho=%1.4f",N,K,FinalTime,x(r,c));
 title(t)            