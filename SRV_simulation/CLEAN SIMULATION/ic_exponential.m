function u0 = ic_exponential(x, N0, A, thickness);
% Exponential decay of initial carrier concentration as a function of the
% absorptance at the excitation wavelength
% thickness is in nm, gets converted to cm to be consistent with units of
% SRV
%   parameter alpha is a function of excitation wavelength
alpha = (A / ((log10((exp(1)))* (thickness*1e-7)))); % A needs to be changed for exciation wavelength, x is the thickness
nc_init_func = N0 * exp(-alpha * x);

u0 = nc_init_func;
end

