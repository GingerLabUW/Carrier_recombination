function [k_nr_fit, k_b_fit] = fit_rate_equation(k_nr_guess, k_b_guess, t, n_meas)
% Fit the rate equation dn/dt = -k_nr * n - k_b * (n^2) to experimental data
% using the guessed values for k_nr and k_b as the starting point for the optimization

% Define the rate equation as a function handle
rate_eqn = @(t, n, k_nr, k_b) -k_nr * n - k_b * (n^2);

% Use the guessed values for k_nr and k_b as the starting point for the optimization
params0 = [k_nr_guess, k_b_guess];

% Define an anonymous function that calculates the sum of squares of the errors
% between the modeled n and the measured n_meas, given values for k_nr and k_b
err_fun = @(params) sum((ode45(@(t,n) rate_eqn(t, n, params(1), params(2)), t, n_meas(1)) - n_meas).^2);

% Use lsqnonlin to fit the rate equation to the experimental data
options = optimoptions('lsqnonlin', 'Display', 'iter');
params_fit = lsqnonlin(err_fun, params0, [], [], options);

% Extract the fitted values for k_nr and k_b
k_nr_fit = params_fit(1);
k_b_fit = params_fit(2);

% Solve the rate equation using the fitted values for k_nr and k_b
n_fit = ode45(@(t,n) rate_eqn(t, n, k_nr_fit, k_b_fit), t, n_meas(1));

end

