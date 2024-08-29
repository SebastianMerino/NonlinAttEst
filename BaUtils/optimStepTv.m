function [A,B,C,D] = optimStepTv(Ja,Jb,DP,rk,alphaCk,betaCk,mu1,mu2,m,n,tol,mask)
% Solver for SLD with Isotropic Total Variation regularization for ACS 
% and Tikhonov regularization for BSC

jcb = [Ja,Jb];
[u,~] = cgs(jcb'*jcb,-jcb'*rk);
A = u(1:m*n);
B = u(m*n+1:end);
C = DP*(alphaCk + A);
D = DP*(betaCk + B);

fid(1) = 1/2*(norm( Ja*A + Ja*B + rk ))^2;
reg(1) = mu1*TVcalc_isotropic(C,m,n,mask) + ...
    mu2*TVcalc_isotropic(D,m,n,mask);
comp(1) = 1/2*norm(C-DP*(alphaCk+A));
ite  = 0;
error = 1;

while abs(error) > tol && ite <= 40
    ite = ite + 1;

    lambda = 1;
    A = solve2norm2(Ja,-(Jb*B+rk),-DP,-(-DP*alphaCk+C), lambda, tol,A);
    B = solve2norm2(Jb,-(Ja*A+rk),-DP,-(-DP*betaCk+D), lambda, tol,B);
    [C,costC] = IRLS_TV(DP*alphaCk+DP*A,speye(m*n),mu1/lambda,m,n,tol,[],mask);
    [D,costD] = IRLS_TV(DP*betaCk+DP*B,speye(m*n),mu2/lambda,m,n,tol,[],mask);
    
    fid(ite+1) = 1/2*(norm( Ja*A + Ja*B + rk ))^2;
    reg(ite+1) = mu1*TVcalc_isotropic(C,m,n,mask) + ...
    mu2*TVcalc_isotropic(D,m,n,mask);
    comp(ite+1) = 1/2*norm(C-DP*(alphaCk+A))^2*lambda;

    if mod(ite,5)==1
        estACtv = reshape(C/5^2,[m,n]);
        estBAtv = reshape(2*(D-1),[m,n]);

        figure('Units','centimeters', 'Position',[5 5 20 10])
        font = 9;
        tiledlayout(1,2)
        nexttile,
        imagesc((1:n)*0.6047,(1:m)*0.6120,estACtv); colorbar;
        clim([0 0.2]);
        title('\alpha_0 in \alpha(f) = \alpha_0 \times f^2 dB/cm')
        axis image
        colormap turbo; colorbar;
        set(gca,'fontsize',font)
        xlabel('Lateral distance (mm)');
        ylabel('Depth (mm)');
        set(gca,'FontSize',font);

        nexttile; imagesc((1:n)*0.6047,(1:m)*0.6120,estBAtv); colorbar;
        clim([5 10]);
        title('B/A');
        axis image
        colormap turbo; colorbar;
        set(gca,'fontsize',font)
        xlabel('Lateral distance (mm)');
        ylabel('Depth (mm)');
        set(gca,'FontSize',font);
        pause(0.1)

        continue
    end
    error = abs(fid(ite+1) + reg(ite+1) - fid(ite) - reg(ite));

end
% figure,tiledlayout(3,1)
% nexttile, plot(fid),
% nexttile, plot(reg)
% nexttile, plot(comp),
end




function u = solve2norm2(A,b,C,d, lambda, tol, u0)
newMat = A'*A + lambda*(C'*C);
newVec = A'*b + lambda*C'*d;
[u,~] = pcg(newMat,newVec,tol,200,[],[],u0);
end
