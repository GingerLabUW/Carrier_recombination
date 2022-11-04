
t = [0 3 5 8 9.5 11.5 14 16 18 20 25 27];
x = [0.0904 0.1503 0.2407 0.3864 0.5201 0.6667 0.8159 0.9979 1.0673 1.1224 1.1512 1.2093]'; 
s = [9.0115 8.8088 7.9229 7.2668 5.3347 4.911 3.5354 1.4041 0 0 0 0]';
p = [0.0151 0.0328 0.0621 0.1259 0.2949 0.3715 0.4199 0.522 0.5345 0.6081 0.07662 0.7869]';
figure
plot(t,[x,s,p],'o-')

xt0 = [0.1;9;0.1]; % read off from data plot [x0 ,S0, p0]';   
p0 = [umax; ks; Yxs; Yps; xt0]; 
Ypred = paramfun1(p0,t); % measurement predictions

pfit = lsqcurvefit(@paramfun1,p0,t,[x,s,p])

umax = pfit(1);  
ks = pfit(2);
Yxs = pfit(3); 
Yps = pfit(4);
f = @(t,a) [umax*a(1)*a(2)/(ks + a(2)); -(1/Yxs)*umax*a(1)*a(2)/(ks + a(2)); (1/Yps)*umax*a(1)*a(2)/(ks + a(2))];
xt0 = pfit(5:7); 
[tspan,a] = ode45(f,[0 25],xt0);    
figure
plot(tspan,a(:,1),tspan,a(:,2),tspan,a(:,3), '-', t,x,'o',t,s,'+',t,p,'s')

function ypred = paramfun1(p,t);
umax = p(1); 
ks = p(2); 
Yxs = p(3); 
Yps = p(4); % use disperse here 
xt0 = p(5:7); % initial conditions 
f = @(t,a) [umax*a(1)*a(2)/(ks + a(2)); -(1/Yxs)*umax*a(1)*a(2)/(ks + a(2)); (1/Yps)*umax*a(1)*a(2)/(ks + a(2))];
[~,ypred] = ode45(f,t,xt0);
end 