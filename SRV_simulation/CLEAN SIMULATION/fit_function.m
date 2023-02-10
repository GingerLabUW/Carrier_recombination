% Define the fitting function
function err = fit_function(params, t_meas, y_meas)
    k_nr = params(1);
    k_b = params(2);
    
    % Solve the rate equation using ode45
    [~, n] = ode45(@(t, n) rate_equation(t, n, k_nr, k_b), t_meas, y_meas(1));
    y_fit = n;
    
    % Calculate the error between the fit and the experimental data
    err = sum((y_fit - y_meas).^2);
    fprintf("err = %f\n", err);
end