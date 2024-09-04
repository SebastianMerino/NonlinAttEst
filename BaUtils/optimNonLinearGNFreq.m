function u = optimNonLinearGNFreq(z,f,bzf, C,d,rho, tol,u0)
%   Optimizes F(u) = 1/2 * ||A(u) - b||^2 + rho/2 * ||C*u - d||^2
%   where A(u) is a non-linear function characterized by the functions
%   model.m and jacobian.m, using the Gauss-Newton approach
% 
%   By Sebastian Merino on August 23rd, 2024

u = u0;

% First residual
res = modelFreq(u,z,f) - bzf(:);
loss = [];
loss(1) = 0.5*norm(res)^2 + rho/2 * norm(C*u - d)^2;

ite = 1;
maxIte = 200;
while true
    jcb = jacobianFreq(u,z,f);
    [step,~] = pcg(jcb'*jcb + rho*(C'*C), -jcb'*res -rho*C'*(C*u-d),tol,200);
    u = u + step;

    % if mod(ite,10)==1
    %     m = 72; n = 63;
    %     alphaArr = reshape(u(1:m*n),[m,n]);
    %     betaArr = reshape(u(1+m*n:end),[m,n]);
    % 
    %     estACtv = alphaArr/5^2/100*8.686;
    %     estBAtv = 2*(betaArr-1);
    % 
    %     fprintf('AC: %.3f +/- %.3f\n', mean(estACtv(:), 'omitnan'), ...
    %     std(estACtv(:), [] ,'omitnan'));
    %     fprintf('B/A: %.2f +/- %.2f\n', mean(estBAtv(:), 'omitnan'), ...
    %     std(estBAtv(:), [] ,'omitnan'));
    % 
    %     figure('Units','centimeters', 'Position',[5 5 20 10])
    %     font = 9;
    %     tiledlayout(1,2)
    %     nexttile,
    %     imagesc((1:n)*0.6047,(1:m)*0.6120,estACtv); colorbar;
    %     clim([0 0.2]);
    %     title('\alpha_0 in \alpha(f) = \alpha_0 \times f^2 dB/cm')
    %     axis image
    %     colormap turbo; colorbar;
    %     set(gca,'fontsize',font)
    %     xlabel('Lateral distance (mm)');
    %     ylabel('Depth (mm)');
    %     set(gca,'FontSize',font);
    % 
    %     nexttile; imagesc((1:n)*0.6047,(1:m)*0.6120,estBAtv); colorbar;
    %     clim([5 10]);
    %     title('B/A');
    %     axis image
    %     colormap pink; colorbar;
    %     set(gca,'fontsize',font)
    %     xlabel('Lateral distance (mm)');
    %     ylabel('Depth (mm)');
    %     set(gca,'FontSize',font);
    %     pause(0.1)
    % end

    res = modelFreq(u,z,f) - bzf(:);
    loss(ite+1) = 0.5*norm(res)^2 + rho/2 * norm(C*u - d)^2;
    if abs(loss(ite+1) - loss(ite))<tol || ite == maxIte, break; end
    % if abs(norm(step))<tol || ite == maxIte, break; end
    ite = ite + 1;
end

% figure,plot(loss)


end