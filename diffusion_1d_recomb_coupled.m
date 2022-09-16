
kdeeptrap = 1e-8;
kmono = 1e6;
kbimol = 7.8e-11;
%diff_coeff = 0.042;
d_array = 0:0.02:0.1;

pumpfwhm = 500e-7;
ncarrierinit = 2.5e12;
ndeeptrap_init = 1e15;



x = linspace(0,5e-4,100);
xsize = size(x);
xcent = x(xsize(2))/2;

xstartidx = find_nearest(x,xcent-(pumpfwhm));
xstopidx = find_nearest(x,xcent+(pumpfwhm));

x_use = x(xstartidx:xstopidx);

t = 0:1e-9:1e-6;
tsize = size(t);
d_array_size = size(d_array);
pl_array = zeros(d_array_size(2),tsize(2));

tmatx = zeros(d_array_size(2),tsize(2));



for i = 1:d_array_size(2)
    tmatx(i,:)=t;
end

parfor i = 1:d_array_size(2)
    m = 0;

    diff_coeff = d_array(i);
    sol = pdepe(m,@(x,t,u,DuDx)diff1dpde(x,t,u,DuDx, kdeeptrap,kmono,kbimol,diff_coeff),@(x)diff1dic(x, pumpfwhm, xcent, ncarrierinit, ndeeptrap_init),@diff1dbc,x,t);
    %%solnodiff = pdepe(m,@nodiff1dpde,@diff1dic,@diff1dbc,x,t);

    u1 = sol(:,:,1);
    %%u1nodiff = solnodiff(:,:,1);



    u1_use = u1(:,xstartidx:xstopidx);
    %%u1_nodiffuse =  u1nodiff(:,xstartidx:xstopidx);
%     figure
%     surf(x./(1e-4),t ./ (1e-9),u1) %'EdgeColor','none'
%     title('N(x,t)')
%     xlabel('Distance (um)')
%     ylabel('Time (ns)')
%     zlabel('N(x,t) (cm-1)')
%     colorbar;
    % 
    % figure
    % surf(x_use./(1e-4),t ./ (1e-9),u1_use) %'EdgeColor','none'
    % title('N(x,t)')
    % xlabel('Distance (um)')
    % ylabel('Time (ns)')
    % zlabel('N(x,t) (cm-1)')
    % colorbar;


    n_integrate = zeros(1,tsize(2));
    %n_integrate_nodiff = zeros(1,tsize(2));
    for j = 1:tsize(2);
        n_int = trapz(x_use, u1_use(j,:));
        %%n_int_nodiff = trapz(x_use, u1_nodiffuse(i,:));

        n_integrate(1,j) = n_int;
        %%n_integrate_nodiff(1,i) = n_int_nodiff;
    end

    pltest = plfun(t, n_integrate, kmono, kbimol);
    %%pltest_nodiff = plfun(t, n_integrate_nodiff);

    pltest_norm = pltest / (pltest(1));
    %%pltest_nodiffnorm = pltest / (pltest_nodiff(1));
    pl_array(i,:) = pltest_norm;
end
    

figure
xlabel('Time (ns)', 'FontSize', 16)
ylabel('PL (Norm.)', 'FontSize', 16)
set(gca,'linewidth',3, 'fontweight', 'bold')

hold on
for i = 1:d_array_size(2)
    pl_array_plot = pl_array(i,:);
    plot(t./(1e-9), pl_array_plot, 'LineWidth', 5);

    
end

% d_array_strings = zeros(1,d_array_size(2));
% for i=1:d_array_size(2)
%     d_array_string = num2str(d_array(1,i));
%     d_array_strings(1,i)=d_array_string
% end
d_array_strings2 = num2str(d_array');

figure
xlabel('Time (ns)', 'FontSize', 16)
ylabel('PL (Norm.)', 'FontSize', 16)
set(gca,'linewidth',3, 'fontweight', 'bold')
hold on
plot(tmatx'./(1e-9), pl_array', 'LineWidth', 5);
legend(num2str(d_array_strings2))


