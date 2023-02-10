
function [c,f,s] = mono_recomb_pde(x,t,u,DuDx, kr, kb, d)
c = 1; 
f = d*DuDx; 
s = -kr*u - kb*(u^2);





% function [c,f,s] = diff1dpde(x,t,u,DuDx, kdt, km, kb, d)
% %kdt is the trapping rate of "deep traps" (dt = deep trap)
% %kdt = 1e-8;
% %km = 1e6;
% %kb = 7.8e-11;
% %d = 0.042;
% c = [1; 1]; 
% f = [d; 0] .* DuDx; 
% s1 = (-kdt * u(2)*u(1))-(km*u(1))-(kb*(u(1)^2));
% s2 = -kdt * u(2)*u(1);
% s = [s1; s2];  