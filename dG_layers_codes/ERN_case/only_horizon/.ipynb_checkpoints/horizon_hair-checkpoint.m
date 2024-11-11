psi_tail = psi_extr(init1:end);
tau_paper = tarr(init1:end) + xL;

hh = 0.5.*tau_paper.*psi_tail;
plot(tau_paper,hh)
saveas(gcf,'horizon_hair.png')