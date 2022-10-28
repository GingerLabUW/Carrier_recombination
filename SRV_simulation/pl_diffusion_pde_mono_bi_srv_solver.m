function [pl_x_integrated] = pl_diffusion_pde_mono_bi_srv_solver(x,t,u,DuDx, kr, kb, d)


sol = pdepe(0,@(x,t,u,DuDx)mono_bi_recomb_pde(x,t,u,DuDx, mono_recomb_coeff,bi_recomb_coeff,diff_coeff),@(x)ic_exponential(x, init_carrier_density, absorbance, film_thickness),@(xl,ul,xr,ur,t)SRV_fixed_bc(xl,ul,xr,ur,t, srv_1, srv_2),x,t_simulate);
% First input correspond to m values it defines symmetry of the parabolic elliptic pde for this equation is always equal 0 (https://www.mathworks.com/help/matlab/ref/pdepe.html#mw_077c5e49-ee92-4c2b-a988-c17d8d564362) section pdefun

%Below is for loop to integrate solution over x for each time point
    
for j = 1:t_size(2); 
    n_int = trapz(x, sol(j,:)); %integrate over x for time point j
    n_integrate(i,j) = n_int; %store in array        
end

%Calculate spatially integrated PL decay from solution to time dependent carrier density

pl_x_integrated = pl_function_mono_bi(t_simulate,n_integrate(i,:), mono_recomb_coeff, bi_recomb_coeff);


