function [u,G] = IRLS_TV_Dy(b,A,mu,M,N,tol,mask,DP)
% Optimizes the following cost function :
%   0.5*||A*u(:)-b||_2^2 + mu*TV(DP*u)
% Inputs: 
%       b               vector containing measurements
%       A               matrix for linear system of eq
%       mu              regularization parameter
%       M,N             image size of u
%       tol             tolerance for error
%  
% Outputs:
%       u               vector of image samples, size MN
%       G               vector containing cost function for each iteration
%
% Modified from a function made by A. Coila
% Matrix DP and docs added by Sebastian Merino

AtA = A'*A;
Atb = A'*b;

% P = sparse(triu(ones(M))')./(1:M)';
% P = kron(speye(N),P);
% Pinv = inv(P);

D = spdiags([-ones(M,1) ones(M,1)], [0 1], M,M+1);
D(:,end) = [];
D(M,M) = 0;
Dy = kron(speye(N),D);

D = spdiags([-ones(N,1) ones(N,1)], [0 1], N,N+1);
D(:,end) = [];
D(N,N) = 0;
Dx = kron(D,speye(M));

Dx = Dx*DP;
Dy = Dy*DP;
D = [Dx' Dy']';

ite = 0;
error = 1;

%[u,~] = cgs(AtA+mu*(D')*D,Atb);
[u,~] = pcg(AtA, Atb);
G(1) = 1/2*(norm( (b - A*u) ))^2 + mu*TVcalc_isotropic(DP*u,M,N,mask);
while error > tol && ite < 200
    
    ite = ite + 1;
    Dh = Dx*u;
    Dv = Dy*u;
    
    vksquare = Dh.^2 + Dv.^2;
    vksquare = vksquare(:);
    
    eps = 0.1; % 0.03
    P = sqrt(vksquare + eps^2);
    P = 1./P;
    
    P = P(:);
    omega = speye(M*N);   % sparse diagonal of just ones
    omega = spdiags(P,0,omega);   % sparse diagonal of P values instead of ones.
    W = kron(speye(2),omega);
    
    [u,~] = pcg(AtA + mu*D'*W*D, Atb,tol,200,[],[],u);
    % [u,~] = cgs(AtA + mu*D'*W*D, Atb,tol,200,[],[],u);
    G(ite+1,1) = 1/2*(norm( (b - A*u) ))^2 + mu*TVcalc_isotropic(DP*u,M,N,mask);
    error = abs(G(ite+1) - G(ite));
    
end

% figure(909); plot(1:length(G),G);

end
