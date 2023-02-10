% Define the function to fit the model to the experimental data
function err = fit_pdepe(params, t_meas, y_meas)

mono_recomb_coeff = params(1);
bi_recomb_coeff = params(2);
diff_coeff = params(3);
init_carrier_density = params(4);
absorbance = params(5);
film_thickness = params(6);
srv_1 = params(7);
srv_2 = params(8);

x = linspace(0, film_thickness, 100); % discretize x
%t = t_meas; % use the same time points as the experimental data

sol = pdepe(0,@(x,t,u,DuDx)mono_bi_recomb_pde(x,t,u,DuDx, mono_recomb_coeff,bi_recomb_coeff,diff_coeff),@(x)ic_exponential(x, init_carrier_density, absorbance, film_thickness),@(xl,ul,xr,ur,t)SRV_fixed_bc(xl,ul,xr,ur,t, srv_1, srv_2), x, t_meas);

y_fit = trapz(x,sol, 2); % integrate over x for each time point

err = sum((y_fit - y_meas).^2); % squared error between the experimental and the fitting data
diff = y_fit(1,1) - y_meas(1,1);
fprintf("err = %f\n", err);
fprintf("diff = %f\n", diff);
end