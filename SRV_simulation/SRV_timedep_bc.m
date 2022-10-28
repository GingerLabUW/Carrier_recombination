% These are the boundary conditions : look here
% (https://www.mathworks.com/help/matlab/ref/pdepe.html#mw_077c5e49-ee92-4c2b-a988-c17d8d564362)
% under bcfunc to learn more

%p is the variable that depends on x,t,u --> in our case that's the SRV 
%q is the diffusion coeffient which multiplies for f whioch is the flux

% we are bascially using this function  for boundary condition: 
% -Sn(x) + D(dn/dx)=0 
% for the top surface : -S1 *ul + D(dn/dx)=0 
% bottom surface : +S2*ur + D(dn/dx)=0
% with: ul,ux = n(x) so the charge carrier concentration in space 
% ql and qr = D(dn/dx)

function [pl,ql,pr,qr] = diff1dbc(xl,ul,xr,ur,t, s1, s2)

pl = -s1*ul; 
ql = 1;     
pr = s2*ur; 
qr = 1;




%Old function
%function [pl,ql,pr,qr] = diff1dbc(xl,ul,xr,ur,t)




%pl = [0; ul(2)]; 
%ql = [1; 0];     
%pr = [0; ur(2)]; 
%qr = [1; 0];    