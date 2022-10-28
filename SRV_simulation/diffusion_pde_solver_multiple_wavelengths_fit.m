
% This is the code to solve the diffusion with SRV boundary conditions 
% First we assign the variables and we refer to them for every function
% between parenteses
%the @ indicates parameters that are only used in the PDE solver
% the name of the variables assigned before the solution do not need to
% match the variables of each function

%remember to add a ; after every line to avoid getting a result cell for
%every code line

%----------------------------------------------------------------

%Here we are reading an ABSORBANCE SPECTRA in .csv to extract the value of
%absorbance at a certain wavelength. Any csv file with first column wvl and
%second column abs values will work

absorbance_file_path = 'C:\Users\Margherita\OneDrive - UW\Documents\DATA\EDA_additive_work\09_17_21_UV_Vis_stab\1day\BR25.csv';

absorbance_data = readtable(absorbance_file_path,'NumHeaderLines',1); %change the number after 'NumHeaderLines' to account for the rows without usable data in your csv file

wavelengths = [400,650,730]; %excitation wavelength in nm

%---------------------------------------------------------------------------

%Here we put the TIME RESOLVED PHOTOLUMINESCENCE DATA
tcspc_data_file_path = 'simulated_wavelength_dependent_pl.csv'; % In this case this is the simulated PL data that is already in the folder at each wavelength selected from the absorbance spectra

tcspc_data = table2array(readtable(tcspc_data_file_path)); % use table2array command to make a data array otherwise it will read it as a table (similar to pandas list Vs np array basically)

%-----------------------------------------------------------------------------------
%here we put the list of the INITIAL VALUES OF THE VARIABLES
%constant values
diff_coeff = 0.004; %units of cm^2/s
mono_recomb_coeff = 1e6; %units of s^-1
bi_recomb_coeff = 4e-11; %units of cm/s
init_carrier_density = 1e13; %units of cm^-1
srv_1 = 1000; %units of cm/s 
srv_2 = 1000; %units of cm/s 
film_thickness = 1000; %in nm

% distance array: x values are a linear sequence (= linspace in python)
xmesh_step = 1; % distance step in nm. It is going to be in the middle of the min and max of the array (min:step:max) shpwn below
x = 0:(xmesh_step*1e-7):(film_thickness*1e-7); %units of cm, the lower and upper limit are the first and last value and the range is in the middle

%time array
t_step = 1e-10; % time step in s
t_simulate = 0:t_step:1e-7; %linspace(1*(10^-6), 5*(10^-6), 501) %units of s

%--------------------------------------------------------------------------------------------------------------
%Now we are goin to do a FOR LOOP to SIMULATE DIFFUSION at each
%wvalength of exciatiton

%here we define the function for the size of the excitation wavlength array
%and time array. We are going to use these tyo iterate in the for loop
wavelengths_size = size(wavelengths); %Size of wavelength array for for loop iteration
t_size = size(t_simulate); % find size of time array

%Here we create the empty datasets that are going to be filled with the
%results of the for loop (like np.zeros in python)

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
    figure
    plot(x, ic_exponential(x, init_carrier_density, absorbance, film_thickness));
    
    %Here is the ACTUAL DIFFUSION PROFILE SIMULATION
    sol = pdepe(0,@(x,t,u,DuDx)mono_bi_recomb_pde(x,t,u,DuDx, mono_recomb_coeff,bi_recomb_coeff,diff_coeff),@(x)ic_exponential(x, init_carrier_density, absorbance, film_thickness),@(xl,ul,xr,ur,t)SRV_fixed_bc(xl,ul,xr,ur,t, srv_1, srv_2),x,t_simulate);
    % First input correspond to m values it defines symmetry of the parabolic elliptic pde for this equation is always equal 0 (https://www.mathworks.com/help/matlab/ref/pdepe.html#mw_077c5e49-ee92-4c2b-a988-c17d8d564362) section pdefun
    %the @ indicates parameters that are only used in the PDE solver

    %OK so what we have now is a distribution of the carrier density in
    %time and space but to look at it in 2D we can integrate the values of
    %carrier density in x for each time point

    %Below is for loop to integrate solution over x for each time point
        
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

% 
% figure
% mesh(x, t_simulate, sol) % to get 3D graphic

figure
set(gca, 'YScale', 'log') %semilog y scale
hold on % this is used to show all the figures at once
for i = 1:wavelengths_size(2);
%hold all
%semilogy(t_simulate, pl_x_integrated);
    plot(t_simulate, pl_norm_integrate(i,:));
end
% 

% simulated_data = zeros(t_size(2), wavelengths_size(2)+1);
% simulated_data(:,1) = t_simulate;
% simulated_data(:,2:end) = pl_integrate.';
% writematrix(simulated_data,'simulated_wavelength_dependent_pl.csv') 


% % 
% % figure
% % plot(x, sol(1,:), 'LineWidth', 5); % to get a single line profile of the
% % change in carrier density in space
% % 
% % figure
% % plot(x, sol(10,:), 'LineWidth', 5);
% % 
% % figure
% % plot(x, sol(1000,:), 'LineWidth', 5);
% 
% 
