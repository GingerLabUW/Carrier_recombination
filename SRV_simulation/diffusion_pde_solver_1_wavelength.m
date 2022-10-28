
% This is the code to solve the diffusion with SRV boundary conditions 
% First we assign the variables and we refer to them for every function
% between parenteses
%the @ indicates parameters that are only used in the PDE solver
% the name of the variables assigned before the solution do not need to
% match the variables of each function

%remember to add a ; after every line to avoid getting a result cell for
%every code line

absorbance_file_path = 'C:\Users\Margherita\OneDrive - UW\Documents\DATA\EDA_additive_work\09_17_21_UV_Vis_stab\1day\BR25.csv';

absorbance_data = readtable(absorbance_file_path,'NumHeaderLines',1);

wavelengths = 730; %excitation wavelength in nm

absorbance = interp1(table2array(absorbance_data(:, 1)),table2array(absorbance_data(:, 2)),wavelength); % value of absorbance at the wavlength of excitation
% 
diff_coeff = 0.04; %units of cm^2/s
mono_recomb_coeff = 1e6; %units of s^-1
bi_recomb_coeff = 4e-10; %units of cm/s
init_carrier_density = 1e14; %units of cm^-1
srv_1 = 1000; %units of cm/s 
srv_2 = 1000; %units of cm/s 

film_thickness = 1000; %in nm

xmesh_step = 1; %in nm
x = 0:(xmesh_step*1e-7):(film_thickness*1e-7); %units of cm, the lower and upper limit are the first and last value and the range is in the middle


t_step = 1e-10;
t_simulate = 0:t_step:1e-6; %linspace(1*(10^-6), 5*(10^-6), 501) %units of s


%Plot initial condition for excitation profile
figure
plot(x, ic_exponential(x, init_carrier_density, absorbance, film_thickness));

sol = pdepe(0,@(x,t,u,DuDx)mono_bi_recomb_pde(x,t,u,DuDx, mono_recomb_coeff,bi_recomb_coeff,diff_coeff),@(x)ic_exponential(x, init_carrier_density, absorbance, film_thickness),@(xl,ul,xr,ur,t)SRV_fixed_bc(xl,ul,xr,ur,t, srv_1, srv_2),x,t_simulate);
% First input correspond to m values it defines symmetry of the parabolic elliptic pde for this equation is always equal 0 (https://www.mathworks.com/help/matlab/ref/pdepe.html#mw_077c5e49-ee92-4c2b-a988-c17d8d564362) section pdefun

%Below is for loop to integrate solution over x for each time point

t_size = size(t_simulate); % find size of time array
n_integrate = zeros(1,t_size(2)); %create empty array to store x integrated solutions for each time point
    
for j = 1:t_size(2);
    n_int = trapz(x, sol(j,:)); %integrate over x for time point j
    n_integrate(1,j) = n_int; %store in array        
end

figure
plot(t_simulate,n_integrate/max(n_integrate))

%Calculate spatially integrated PL decay from solution to time dependent carrier density

pl_x_integrated = pl_function_mono_bi(t_simulate,n_integrate, mono_recomb_coeff, bi_recomb_coeff);

pl_x_integrated_norm = pl_x_integrated / max(pl_x_integrated);


simulated_noise = sqrt(poissrnd(pl_x_integrated * t_step));


noisy_pl = simulated_noise;



figure
mesh(x, t_simulate, sol) % to get 3D graphic

figure
%hold all
%semilogy(t_simulate, pl_x_integrated);
semilogy(t_simulate, noisy_pl);
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
