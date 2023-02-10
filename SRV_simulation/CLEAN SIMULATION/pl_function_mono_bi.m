function [pl] = pl_function_mono_bi(t,n, km, kb);
pl = (km.*n) + (kb.*(n.^2));