function [B,C] = optimAdmmTvDy(A1,A2,b,mu1,mu2,m,n,tol,mask)
% Solver for SLD with Isotropic Total Variation regularization for ACS 
% and Tikhonov regularization for BSC

A = [A1 A2];
AtA = A'*A;
Atb = A'*b;
[u,~] = cgs(AtA,Atb);
B = reshape(u(1:end/2),m,n);
C = reshape(u(end/2+1:end),m,n);

B = B(:);
C = C(:);
D = 0;
v = 0;

F(1) = 1/2*(norm( b - A1*B - A2*C ))^2 + ...
    mu1*TVcalc_isotropic(B,m,n,mask) + mu2*sum(abs(C(:)),'all');
ite  = 0;
error = 1;

while abs(error) > tol && ite < 20
    ite = ite + 1;
    
    rho = 1;
    % First part of ADMM algorithm: B
    [B,costB] = IRLS_TV_Dy(b-A2*C-D-v,A1,mu1/rho,m,n,tol,mask);
    
    % Second part of ADMM algorithm: C
    [C,costC] = IRLS_TV_Dy(b-A1*B-D-v,A2,mu2/rho,m,n,tol,mask);
    
    % Third part of ADMM algorithm: D
    % Least squares: 1/2*||D||_2^2 + rho/2*||D-w||_2^2
    w = b - A1*B - A2*C - v;
    D = (rho/(rho+1))*w;
    
    % Fourth part of ADMM algorithm: v
    v = v + A1*B + A2*C + D - b;
    F(ite+1,1) = 1/2*(norm( b - A1*B - A2*C ))^2 + ...
        mu1*TVcalc_isotropic(B,m,n,mask) + mu2*TVcalc_isotropic(C,m,n,mask);
    error = abs(F(ite+1) - F(ite));
end

end


