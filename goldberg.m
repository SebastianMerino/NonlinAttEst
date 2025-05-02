f = 6e6;
beta = 1 + 12/2;
P = 400e3;
rho = 1000;
c0 = 1540;

NptodB = db(exp(1));
alpha = 0.08*(f/1e6)^2;
alpha = alpha/NptodB*100; % Np/m
GN = 2*pi*f*P*beta/c0^3/rho/alpha