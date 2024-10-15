function jcb = jacobianFreqPL(u,z,freq)
% Gets the jacobian matrix, of the function specified in modelFreq.m
% given u = [a0' gm' b']', the column vector z and the row vector freq.
%   a0 is in Np/m/MHz^gm,
%   gm is adimensional
%   b is adimensional, 
%   z is in m, 
%   freq is in MHz 
% By Sebastian Merino

m = length(z);     % number of rows
n = length(u)/3/m; % number of cols
p = length(freq);  % number of frequency points

f(1,1,:) = freq;
alpha0Map = reshape(u(1:m*n),m,n);
gammaMap = reshape(u(m*n+1:2*m*n),m,n);
acMap = alpha0Map.*(f.^gammaMap);
betaMap = reshape(u(2*m*n+1:end),m,n);

dAalpha0 = betaMap.*( (2*acMap.*z+1).*exp(-2*acMap.*z) - 1)....
    ./(acMap.^2.*z) .*(f.^gammaMap);
dAgamma = betaMap.*( (2*acMap.*z+1).*exp(-2*acMap.*z) - 1)....
    ./(acMap.^2.*z) .*( alpha0Map .* log(f) .* (f.^gammaMap) );
dAbeta = (1./(acMap.*z)).*(1 - exp(-2*acMap.*z));

% Indices for sparse matrix
I = (1:m*n*p)';
J = repmat((1:m*n)',p,1);

jcb = [sparse(I,J,dAalpha0(:)),sparse(I,J,dAgamma(:)),sparse(I,J,dAbeta(:))];
end
