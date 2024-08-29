function u = optimNonLinearGN_v2(X,Y, C,d,rho, tol,u0)
%   Optimizes F(u) = 1/2 * ||A(u) - b||^2 + rho/2 * ||C*u - d||^2
%   where A(u) is a non-linear function characterized by the functions
%   model.m and jacobian.m, using the Gauss-Newton approach
% 
%   By Sebastian Merino on August 23rd, 2024

u = u0;

% First residual
res = Y - model_v2(u,X); 
loss = [];
loss(1) = 0.5*norm(res)^2 + rho/2 * norm(C*u - d)^2;

ite = 1;
maxIte = 200;
while true
    jcb = jacobian_v2(u,X);
    [step,~] = pcg(jcb'*jcb + rho*(C'*C), -jcb'*res -rho*C'*(C*u-d),tol,200);
    u = u + step;

    res = Y - model_v2(u,X); % Next residual
    loss(ite+1) = 0.5*norm(res)^2 + rho/2 * norm(C*u - d)^2;
    if abs(loss(ite+1) - loss(ite))<tol || ite == maxIte, break; end
    % if abs(norm(step))<tol || ite == maxIte, break; end
    ite = ite + 1;
end

% figure,plot(loss)


end