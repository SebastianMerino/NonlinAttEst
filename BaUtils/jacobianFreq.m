function jcb = jacobianFreq(u,z,freq, varargin)
% Gets the jacobian matrix, of the function specified in modelFreq.m
% given u = [a' b']', the column vector z and the row vector freq.
%   a is in Np/m/MHz^2, 
%   b is adimensional, 
%   z is in m, 
%   freq is in MHz 
% By Sebastian Merino

m = length(z);     % number of rows
n = length(u)/2/m; % number of cols
p = length(freq);  % number of frequency points

if nargin == 3
    f2(1,1,:) = freq.^2;
else
    f2(1,1,:) = freq.^varargin{1};
end
Z = repmat(z,1,n);
f2Z = f2.*Z;

alphaMap = reshape(u(1:m*n),m,n);
betaMap = reshape(u(m*n+1:end),m,n);

dAalpha = betaMap.*( (2*f2Z.*alphaMap+1).*exp(-2*alphaMap.*f2Z) - 1)....
    ./(alphaMap.^2.*f2Z);
dAbeta = (1./(alphaMap.*f2Z)).*(1 - exp(-2*alphaMap.*f2Z));

% Indices for sparse matrix
I = (1:m*n*p)';
J = repmat((1:m*n)',p,1);

jcb = [sparse(I,J,dAalpha(:)) sparse(I,J,dAbeta(:))];
end
