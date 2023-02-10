% Load in the absorbance data from a .csv file
absorbance_file_path = 'C:\Users\Margherita\OneDrive - UW\Documents\DATA\EDA_additive_work\09_17_21_UV_Vis_stab\1day\BR25.csv';
absorbance_data = readtable(absorbance_file_path,'NumHeaderLines',1); % Change the number after 'NumHeaderLines' to account for the number of rows without usable data in the csv file

% Define the wavelengths (in nm) for the excitation
wavelengths = [400,650,730];

% Load in the tcspc data from a .csv file
tcspc_data_file_path = 'simulated_wavelength_dependent_pl.csv';
tcspc_data = table2array(readtable("640nm_LLI.csv"));

% Define the diffusion coefficient, mono-recombination coefficient, bi-recombination coefficient, initial carrier density, SRV 1 and SRV 2, and film thickness
diff_coeff = 0.004; % Units of cm^2/s
mono_recomb_coeff = 1e6; % Units of s^-1
bi_recomb_coeff = 4e-11; % Units of cm/s
init_carrier_density = 1e13; % Units of cm^-1
srv_1 = 1000; % Units of cm/s
srv_2 = 1000; % Units of cm/s
film_thickness = 1000; % In nm

% Define the x-mesh step
xmesh_step = 1; 
x = 0:(xmesh_step*1e-7):(film_thickness*1e-7);

%time array
t_step = 1e-10; % time step in s
t_simulate = 0:t_step:1e-6; %linspace(1*(10^-6), 5*(10^-6), 501) %units of s

%--------------------------------------------------------------------------------------------------------------
%Now we are goin to do a FOR LOOP to SIMULATE DIFFUSION at each wavelength of exciatiton

%here we define the function for the size of the excitation wavelength array and time array. We are going to use these to iterate in the for loop
wavelengths_size = size(wavelengths); %Size of wavelength array for for loop iteration
t_size = size(t_simulate); % find size of time array

%Here we create the empty datasets that are going to be filled with the results of the for loop (like np.zeros in python)

% **1st array to fill**
n_integrate = zeros(wavelengths_size(2),t_size(2)); %create empty matrix to store x integrated solutions for each time point and wavelength
% **2nd array to fill**
pl_integrate =zeros(wavelengths_size(2),t_size(2));%create empty matrix to store x integrated PL simulations for each time point and wavelength
% **3rd array to fill**
pl_norm_integrate = zeros(wavelengths_size(2),t_size(2)); %create empty matrix to store x integrated normalized PL simulations for each time point and wavelength

%FOR LOOP to create new datasets starts here!

%**1st array to fill**
for i = 1:wavelengths_size(2); %for loop for each excitation wavelength. Remeber we need to add the (2) after size because Matlab creates an array of [1, size] so we need to cale the value of the second index of the array
    wavelength = wavelengths(i); %we call the first wavelength we want to use
    absorbance = interp1(table2array(absorbance_data(:, 1)),table2array(absorbance_data(:, 2)),wavelength); % value of absorbance at the wavelength of excitation. This function is interpolating the absorbance array x and y (first two values) to find the value of absorbance at the wvalength of exciation (third value)

%Plot initial condition for excitation profile
    %figure
    %plot(x, ic_exponential(x, init_carrier_density, absorbance, film_thickness));
    
    %Here is the ACTUAL DIFFUSION PROFILE SIMULATION
    sol = pdepe(0,@(x,t,u,DuDx)mono_bi_recomb_pde(x,t,u,DuDx, mono_recomb_coeff,bi_recomb_coeff,diff_coeff),@(x)ic_exponential(x, init_carrier_density, absorbance, film_thickness),@(xl,ul,xr,ur,t)SRV_fixed_bc(xl,ul,xr,ur,t, srv_1, srv_2),x,t_simulate);
    % First input correspond to m values it defines symmetry of the parabolic elliptic pde for this equation is always equal 0 (https://www.mathworks.com/help/matlab/ref/pdepe.html#mw_077c5e49-ee92-4c2b-a988-c17d8d564362) section pdefun the @ indicates parameters that are only used in the PDE solver

%the next loop integrates the solution of the simulation over x  
    for j = 1:t_size(2); %the number (2) in time size just indicate that we are selecting the size in the [1, size] matrix created by Matlab
        n_int = trapz(x, sol(j,:)); %integrate over x for time point j. The function trapz(x,y) integrates Y with respect to the coordinates or scalar spacing specified by X. We want a trapezoid because it is the closest shape that we get when we slice the exponential decay in time slices
        n_integrate(i,j) = n_int; %FIRST empty array filled here!  
    end

    %**2nd array to fill**

    %Calculate spatially integrated PL decay from solution to time dependent carrier density
    
    pl_x_integrated = pl_function_mono_bi(t_simulate,n_integrate(i,:), mono_recomb_coeff, bi_recomb_coeff);
    pl_integrate(i,:) = sqrt(poissrnd(pl_x_integrated * t_step)) ; %this function is used to reproduce the noise of the actual TRPL decay that increase when the intensity of signal is low
    
    %**3rd array to fill**
    pl_norm_integrate(i,:) = pl_integrate(i,:) / max(pl_integrate(i,:)); % this it to get normalized PL decays

end

%This one is to plot the 3D graph of solved system of differential
%equations

%figure
%mesh(x, t_simulate, sol) % to get 3D graphic

figure
hold on
set(gca, 'YScale', 'log') %semilog y scale
for i = 1:wavelengths_size(2);
hold all
%plot(t_simulate, pl_x_integrated);
%plot(t_simulate, n_integrate(i,:)/max(n_integrate(i,:)));
%plot(t_simulate, pl_norm_integrate(i,:));
end


%---------------------------------------FITTING------------------------------------------

%Load the experimental data from the tcspc_data matrix
% The first column of tcspc_data is the time data, and it is converted to seconds
t_meas = tcspc_data(:,1)*1e-9; % Load the x_meas data
%t_meas = t_simulate;

% The second column of tcspc_data is the intensity data, and it is normalized by dividing by the maximum value
y_meas = tcspc_data(:,2); % Load the y_meas data
%y_meas =  n_integrate(2,:)/max(n_integrate(2,:));

%plot(t_meas, y_meas, '-', 'MarkerSize', 10, 'LineWidth', 2); 

% Set the initial guess for the parameters to be optimized
params0 = [ 3.5392e6 , 4e-9, 0.001, 100000, 2.3, 500, 10000, 100];
lb = [ 3.5392e6, 4e-9, 0.001, 100000, 2.3, 500, 10000, 100];
ub = [ 3.5392e6, 4e-9, 0.001, 200000, 2.3, 500, 10000, 100];

% Set the options for the non-linear least squares optimization algorithm
options = optimoptions('lsqnonlin', 'Display', 'iter',  'TolFun', 1e-40, 'MaxIterations', 1);

% Use lsqnonlin to optimize the parameters to fit the experimental data
params_fit = lsqnonlin(@(params) fit_pdepe(params, t_meas, y_meas), params0, lb, ub, options);

% Solve the partial differential equation using the optimized parameters
sol = pdepe(0,@(x,t,u,DuDx)mono_bi_recomb_pde(x,t,u,DuDx, params_fit(1),params_fit(2),params_fit(3)),@(x)ic_exponential(x,  params_fit(4),params_fit(5),params_fit(6)),@(xl,ul,xr,ur,t)SRV_fixed_bc(xl,ul,xr,ur,t, params_fit(7), params_fit(8)), x, t_meas);

% Integrate the solution along the x dimension to get the total carrier concentration
y_fit = trapz(x,sol,2);


% Set the initial guess for the parameters to be optimized
params0_2 = [ 3.5392e6 , 4e-9, 0.001, 200000, 4, 500, 10000, 100];
lb_2 = [ 3.5392e6, 4e-9, 0.001, 100000, 4, 500, 10000, 100];
ub_2 = [ 3.5392e6, 4e-9, 0.001, 400000, 4, 500, 10000, 100];

% Set the options for the non-linear least squares optimization algorithm
options_2 = optimoptions('lsqnonlin', 'Display', 'iter',  'TolFun', 1e-40, 'MaxIterations', 1);

% Use lsqnonlin to optimize the parameters to fit the experimental data
params_fit_2 = lsqnonlin(@(params) fit_pdepe(params, t_meas, y_meas), params0_2, lb_2, ub_2, options_2);

% Solve the partial differential equation using the optimized parameters
sol_2 = pdepe(0,@(x,t,u,DuDx)mono_bi_recomb_pde(x,t,u,DuDx, params_fit_2(1),params_fit_2(2),params_fit_2(3)),@(x)ic_exponential(x,  params_fit_2(4),params_fit_2(5),params_fit_2(6)),@(xl,ul,xr,ur,t)SRV_fixed_bc(xl,ul,xr,ur,t, params_fit_2(7), params_fit_2(8)), x, t_meas);

% Integrate the solution along the x dimension to get the total carrier concentration
y_fit_2 = trapz(x,sol_2,2);

% Plot the fitting curve on top of the experimental data
hold on
plot(t_meas, y_meas, 'r-', 'MarkerSize', 10, 'LineWidth', 2); 
plot(t_meas, y_fit, '-', 'LineWidth', 2);
plot(t_meas, y_fit_2, '-', 'LineWidth', 2);
xlabel('Time (s)');
set(gca, 'YScale', 'log') %semilog y scale
ylabel('Counts');
legend({'Experimental Data', 'Fitting Curve A =2.3', 'Fitting Curve A =4'});
title('Fitting Curve and Experimental Data');