function mdl = modelFreqLinear(u,z,freq)
% Evaluates the model given u = [a0' a1' b']', the column vector z, and the row
% vector freq.
%   a0 is in Np/m,
%   a1 is in Np/m/MHz
%   b is adimensional, 
%   z is in m, 
%   freq is in MHz 
% Result is a 3D map (rows, columns, freq) stacked in a m*n*p vector
% By Sebastian Merino

m = length(z);     % number of rows
n = length(u)/3/m; % number of cols

f(1,1,:) = freq;
alpha0Map = reshape(u(1:m*n),m,n);
alpha1Map = reshape(u(m*n+1:2*m*n),m,n);
acMap = alpha0Map + alpha1Map.*f;
betaMap = reshape(u(2*m*n+1:end),m,n);

mdl = (betaMap./(z.*acMap)).*(1 - exp(-2*(acMap.*z))) ;
mdl = mdl(:);
end