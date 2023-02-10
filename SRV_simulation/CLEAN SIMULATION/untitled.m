% Load the experimental data for n as a function of t
tcspc_data = table2array(readtable("Br25_TRPL.csv"));
t = tcspc_data(:,1)*1e-9; % time values
n_meas = tcspc_data(:,2)/max(tcspc_data(:,2)); % concentration values

% Define the rate equation function
rate_equation = @(t,n,knr,kr) -knr*n - kr*n.^2;

% Guess values for knr and kr
knr_guess = 1e6;
kr_guess = 4e-11;

% Fit the rate equation to the experimental data
params0 = [knr_guess, kr_guess];
options = optimset('MaxFunEvals', 1e5, 'MaxIter', 1e5);
params_fit = lsqnonlin(@(params) n_meas - deval(ode45(@(t,n) rate_equation(t,n,params(1),params(2)), t, n_meas), params0, [], [], options);

% Get the fitted values of knr and kr
knr_fit = params_fit(1);
kr_fit = params_fit(2);

% Plot the experimental data and the fit
figure
plot(t, n_meas, 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
hold on
plot(t, deval(ode45(@(t,n) rate_equation(t,n,knr_fit,kr_fit), t, n_meas), '-r', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Concentration');
legend('Experimental Data', 'Fit');
title('Experimental Data Fit to Recombination Rate Equation');

% Display the fitted values of knr and kr
disp(['Fitted knr value: ', num2str(knr_fit), ' s^-1']);
disp(['Fitted kr value: ', num2str(kr_fit), ' s^-1']);