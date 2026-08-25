function Dt = Dt_1dperiodic(t, D0, a, kappa)
% 
% Calculates time-dependent cumulative diffusion coefficient D(t)
% from a system of periodic permeable barriers 
% according to Eq. S25 from Suppl Info of 
% Novikov et al PNAS 111, 5088 (2014)
% http://www.pnas.org/cgi/doi/10.1073/pnas.1316944111
% period = a, free diffusivity = D0, permeability = kappa
% 
% (c) Dmitry S Novikov dima@alum.mit.edu

zeta = D0/(kappa*a);  
td = a^2/D0; 
tbar = t/td; 

M = 40; beta = arrayfun(@(m) fzero(@(k) tan(k/2)+zeta*k/2, [(2*m-1)*pi+1e-12, 2*m*pi]), (1:M)'); 
Dt = D0*(1/(1+zeta) + sum((1 - exp(-beta.^2 * tbar))./beta.^2./(1 + zeta + (zeta*beta/2).^2), 1) * 2*zeta^2./tbar);

%figure; hold on; grid on; 
%plot(tbar.^(1/2), Dt/D0, 'k.'); 
%plot(tbar.^(1/2), 1 - 8/3/sqrt(pi)*sqrt(tbar), 'b--');
