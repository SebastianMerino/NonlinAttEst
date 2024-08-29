function jcb = jacobian_v2(theta,x)
mn = size(theta,1)/2; % number of points in image
p = length(x)/mn;

alphaArr = theta(1:mn)'/db(exp(1))*100; % dB/cm -> Np/m
betaArr = 1./theta(mn+1:end)';
zBA = reshape(x,[p,mn]);

dMLfda = -(betaArr).*((2*zBA.*alphaArr+1).*...
    exp(-2*alphaArr.*zBA) - 1)./(alphaArr.^2.*zBA) /db(exp(1))*100;
dMLfdb = -(1./zBA./alphaArr).*(1 - exp(-2*alphaArr.*zBA)).*(-betaArr.^2);

% Indices for sparse matrix
I = reshape(1:mn*p,[mn,p]); % Row indices
J = (1:mn).*ones(p,1);

jcb = [sparse(I,J,dMLfda) sparse(I,J,dMLfdb)];
end
