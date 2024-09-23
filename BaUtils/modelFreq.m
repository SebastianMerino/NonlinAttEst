function mdl = modelFreq(u,z,freq)
% Evaluates the model given u = [a' b']', the column vector z, and the row
% vector freq.
%   a is in Np/m/MHz^2, 
%   b is adimensional, 
%   z is in m, 
%   freq is in MHz 
% Result is a 3D map (rows, columns, freq) stacked in a m*n*p vector
% By Sebastian Merino

m = length(z);     % number of rows
n = length(u)/2/m; % number of cols

f2(1,1,:) = freq.^2;
alphaMap = reshape(u(1:m*n),m,n);
betaMap = reshape(u(m*n+1:end),m,n);

mdl = (betaMap./(z.*alphaMap.*f2)).*(1 - exp(-2*(alphaMap.*f2.*z))) ;
mdl = mdl(:);
end