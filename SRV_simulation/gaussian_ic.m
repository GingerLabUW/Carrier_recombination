
%ic stands for initial condition and below we are expressing the initial
%charge carrier concentration as a gaussian profile


function u0 = diff1dic(x, fwhm, gausscent, nc_init);
%this is a function that defines the initial conditions
%fwhm = 500e-7;
sig = fwhm / (2*sqrt(2*log(2)));
%gausscent = 2.5e-4;
%nc_init = 2.5e12;
nc_init_fun = nc_init*(1/(sig*sqrt(2*pi)))*exp(-((x-gausscent)^2)/(2*(sig^2)));

u0 = nc_init_fun;

% 
% function u0 = diff1dic(x, n_init)
% u0 = n_init*dirac(x)

