function [step,uL] = optimAdmmStepTv(jcb,DP,res,uC,mu1,mu2,m,n,tol,mask)
% Solver for SLD with Isotropic Total Variation regularization for ACS 
% and Tikhonov regularization for BSC

[step,~] = cgs(jcb'*jcb,-jcb'*res);
DPC = blkdiag(DP,DP);
uL = DPC*uC;
I = speye(2*m*n);
Ii = I(:,1:m*n);
Id = I(:,1+m*n:end);

fid(1) = 1/2*(norm( jcb*step + res ))^2;
reg(1) = mu1*TVcalc_isotropic(uL(1:m*n),m,n,mask) + ...
mu2*TVcalc_isotropic(uL(1+m*n:end),m,n,mask);
comp(1) = 1/2*norm(uL - DPC*uC - DPC*step)^2;
w = 0;
ite  = 0;
error = 1;
rho = 1;

while abs(error) > tol && ite <= 20
    ite = ite + 1;

    step = solve2norm2(jcb,-res,-DPC,-uL+DPC*uC-w, rho,tol);
    [alphaL,betaL] = AlterOpti_ADMM(Ii,Id,DPC*step+DPC*uC-w,mu1,mu2,...
        m,n,tol,[mask;mask]);
    uL = [alphaL;betaL];
    w = w - DPC*step + uL - DPC*uC;

    fid(ite+1) = 1/2*(norm( jcb*step + res ))^2;
    reg(ite+1) = mu1*TVcalc_isotropic(uL(1:m*n),m,n,mask) + ...
    mu2*TVcalc_isotropic(uL(1+m*n:end),m,n,mask);
    comp(ite+1) = 1/2*norm(uL - DPC*uC - DPC*step)^2;

    if ite == 20
        estACtv = reshape(uL(1:m*n)/5^2,m,n);
        estBAtv = reshape(2*(uL(1+m*n:end)-1),m,n);

        figure('Units','centimeters', 'Position',[5 5 20 10])
        font = 9;
        tiledlayout(1,2)
        t1 = nexttile;
        imagesc((1:n)*0.6047,(1:m)*0.6120,estACtv); colorbar;
        clim([0 0.2]);
        title('\alpha_0 in \alpha(f) = \alpha_0 \times f^2 dB/cm')
        axis image
        colorbar;
        set(gca,'fontsize',font)
        xlabel('Lateral distance (mm)');
        ylabel('Depth (mm)');
        set(gca,'FontSize',font);

        nexttile; imagesc((1:n)*0.6047,(1:m)*0.6120,estBAtv); colorbar;
        clim([5 10]);
        title('B/A');
        axis image
        colormap pink; colorbar;
        set(gca,'fontsize',font)
        xlabel('Lateral distance (mm)');
        ylabel('Depth (mm)');
        set(gca,'FontSize',font);

        colormap(t1,turbo)
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




function u = solve2norm2(A,b,C,d, lambda, tol)
newMat = A'*A + lambda*(C'*C);
newVec = A'*b + lambda*C'*d;
[u,~] = pcg(newMat,newVec,tol,200);
end
