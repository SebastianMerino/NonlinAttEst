%% B/A estimator with time domain approach
% by: Andres Coila

clear; close all; clc;
addpath(genpath(pwd))
% load('test.mat');

%% Getting system
freq = 5;
radiusDisk = (9)*1e-3;
centerDepth = 22.5e-3;

dx = 3.0236e-04;
dz = 9.0000e-06;
m = 6923; n = 128;
x = (0:n-1)*dx; x = x - mean(x);
z = (0:m-1)'*dz;

%% Generating local maps
NptodB = db(exp(1));
[Xmesh,Zmesh] = meshgrid(x,z);
inc = Xmesh.^2 + (Zmesh-centerDepth).^2 <= (radiusDisk).^2;

bzf = [];
freqVec = [4,5,6];
for iFreq = 1:length(freqVec)
    betaL = ones(size(Xmesh))*(1+6/2);
    betaL(inc) = (1+9/2);
    alphaL = ones(size(Xmesh))*0.1;
    freq = freqVec(iFreq);
    
    figure('Units','centimeters', 'Position',[5 5 20 10])
    font = 9;
    tiledlayout(1,2)
    t1 = nexttile;
    imagesc(x,z,alphaL); colorbar;
    clim([0 0.2]);
    title('\alpha_0 in \alpha(f) = \alpha_0 \times f^2 dB/cm')
    axis image
    colorbar;
    set(gca,'fontsize',font)
    xlabel('Lateral distance (cm)');
    ylabel('Depth (cm)');
    set(gca,'FontSize',font);
    nexttile; imagesc(x,z,(betaL-1)*2); colorbar;
    clim([5 10]);
    title('B/A');
    axis image
    colormap pink; colorbar;
    set(gca,'fontsize',font)
    xlabel('Lateral distance (cm)');
    ylabel('Depth (cm)');
    set(gca,'FontSize',font);
    colormap(t1,turbo)
    pause(0.1)
    
    % Generating cumulative maps
    izBlock = round(z./dz);
    betaC = cumsum(betaL)./izBlock;
    alphaC = cumsum(alphaL)./izBlock;
    
    figure('Units','centimeters', 'Position',[5 5 20 10])
    font = 9;
    tiledlayout(1,2)
    t1 = nexttile;
    imagesc(x*1e2,z*1e2,alphaC); colorbar;
    clim([0 0.2]);
    title('\alpha_0 in \alpha(f) = \alpha_0 \times f^2 dB/cm')
    axis image
    colorbar;
    set(gca,'fontsize',font)
    xlabel('Lateral distance (cm)');
    ylabel('Depth (cm)');
    set(gca,'FontSize',font);
    nexttile; imagesc(x*1e2,z*1e2,(betaC-1)*2); colorbar;
    clim([5 10]);
    title('B/A');
    axis image
    colormap pink; colorbar;
    set(gca,'fontsize',font)
    xlabel('Lateral distance (cm)');
    ylabel('Depth (cm)');
    set(gca,'FontSize',font);
    colormap(t1,turbo)
    pause(0.1)
    
    alphaCnp = alphaC/NptodB*100 * freq^2; % dB/cm -> Np/m
    mzaxis = betaC.*(1 - exp(-2*alphaCnp.*Zmesh) )./alphaCnp./Zmesh + ...
        0.0*randn(size(Xmesh));
    %% Getting system
    % Subsampling parameters
    wl = 1540/5/1e6; % Mean central frequency
    blockParams.blockSize = [10 10]*wl; 
    blockParams.overlap = 0.8;
    blockParams.zlim = [0.3; 5.5]/100;
    blockParams.xlim = [-2; 2]/100;
    
    [bz,xP,zP] = getMeanBlock(mzaxis,x,z,blockParams);
        bzf(:,:,iFreq) = bz;
end

%%
figure,
tiledlayout(1,3)
nexttile,
imagesc(xP*100,zP*100,bzf(:,:,1), [2 13])
title("Measurements b(z,f)")
xlabel('Lateral [cm]')
ylabel('Depth [cm]')
axis image
colorbar
colormap parula

nexttile,
imagesc(xP*100,zP*100,bzf(:,:,2), [2 13])
title("Measurements b(z,f)")
xlabel('Lateral [cm]')
ylabel('Depth [cm]')
axis image
colorbar
colormap parula

nexttile,
imagesc(xP*100,zP*100,bzf(:,:,3), [2 13])
title("Measurements b(z,f)")
xlabel('Lateral [cm]')
ylabel('Depth [cm]')
axis image
colorbar
colormap parula
%% Gauss-Newton with LM
[m,n,p] = size(bzf);
tol = 1e-3;
maxIte = 400;
muAlpha = 0;
muBeta = 0;
% muAlpha = 0; muBeta = 0;
betaIni = 1+(10.5)/2;
alphaIni = 0.08/NptodB*100;

u = [alphaIni*ones(n*m,1);betaIni*ones(n*m,1)];
regMatrix = blkdiag(muAlpha*speye(n*m),muBeta*speye(n*m));

loss = [];
loss(1) = 0.5*norm(modelFreq(u,zP,freqVec) - bzf(:))^2;

ite = 1;
tic 
while true
    jcb = jacobianFreq(u,zP,freqVec);
    res = modelFreq(u,zP,freqVec) - bzf(:);
    [step,~] = pcg(jcb'*jcb + regMatrix,jcb'*-res, tol, 200);
    % step = (jcb'*jcb + regMatrix)\(jcb'*-res);
    u = u + step;

    loss(ite+1) = 0.5*norm(modelFreq(u,zP,freqVec) - bzf(:))^2;
    if abs(loss(ite+1) - loss(ite))<tol || ite == maxIte, break; end
    ite = ite + 1;
end
toc

alphaArr = u(1:n*m);
betaArr = u(n*m+1:end);

estAClm = reshape(alphaArr*NptodB/100,[m,n]);
estBAlm = reshape(2*(betaArr-1),[m,n]);

fprintf('AC: %.3f +/- %.3f\n', mean(estAClm(:), 'omitnan'), ...
    std(estAClm(:), [] ,'omitnan'));
fprintf('B/A: %.2f +/- %.2f\n', mean(estBAlm(:), 'omitnan'), ...
    std(estBAlm(:), [] ,'omitnan'));

figure('Units','centimeters', 'Position',[5 5 10 5]), 
plot(loss, 'LineWidth',2)
xlim([2 length(loss)])
xlabel('Number of iterations')
ylabel('Loss')

figure; imagesc(xP*1e3,zP*1e3,estAClm); colorbar;
clim([0 0.2]);
%title('ACS (dB/100/MHz)');
title('\alpha_0 in \alpha(f) = \alpha_0 \times f^2 dB/100')
font = 20;
axis image
colormap turbo; colorbar;
%clim([0 12])
%font = 20;
set(gca,'fontsize',font)
xlabel('Lateral distance (mm)');
ylabel('Depth (mm)');
%title("\bf BoA map");
set(gca,'FontSize',font);
%ylim([5 55])


figure; imagesc(xP*1e3,zP*1e3,estBAlm); colorbar;
clim([5 10]);
title('B/A');
title(['B/A = ',num2str(median(estBAlm(:),'omitnan'),'%.1f')]);
font = 20;
axis image
colormap pink; colorbar;
set(gca,'fontsize',font)
xlabel('Lateral distance (mm)');
ylabel('Depth (mm)');
%title("\bf BoA map");
set(gca,'FontSize',font);
pause(0.5)


%% Inversion, Coila's implementation
zBlock = zP; xBLock = xP; 
muLocal = 0.01;
dzBlock = zBlock(2)-zBlock(1);
izBlock = round(zBlock./dzBlock);

factorq = izBlock(1)./izBlock ;
estBAcum = estBAlm - estBAlm(1,:).*factorq; 
% estBAcum = estBAtv - estBAtv(1,:).*factorq;
estBAcum = estBAcum(2:end,:);


P = sparse(tril(ones(m-1)));
P = P./izBlock(2:end);
P = kron(speye(n),P);
estBAinst = IRLS_TV(estBAcum(:),P,muLocal,m-1,n,tol,[],ones((m-1)*n,1));
estBAinst = reshape(estBAinst,m-1,n);

figure; imagesc(x*1e3,z*1e3,estBAinst); colorbar;
clim([5 10]);
title('B/A');
font = 20;
axis image
colormap pink; colorbar;
set(gca,'fontsize',font)
xlabel('Lateral distance (mm)');
ylabel('Depth (mm)');
hold on
rectangle('Position',[0-radiusDisk,centerDepth-radiusDisk,...
2*radiusDisk,2*radiusDisk]*1000, 'Curvature',1,...
'EdgeColor','b', 'LineStyle','--', 'LineWidth',2)
hold off
set(gca,'FontSize',font);
pause(0.1)


%% ADMM
% Hyperparameters
tol = 1e-3;
muAlpha = 1e-1; muBeta = 1e-1; % USED IN PRESENTATION
% muAlpha = 1e-3; muBeta = 1e-3;
% muAlpha = 10; muBeta = 10;
rho = 0.01;
maxIte = 50;

% Initialization
betaIni = 1+(7.5)/2;
alphaIni = 0.08/NptodB*100;
u = [alphaIni*ones(n*m,1);betaIni*ones(n*m,1)];
v = [alphaIni*ones(n*m,1);betaIni*ones(n*m,1)];
w = zeros(size(u));

% Auxiliar
dzBlock = zBlock(2)-zBlock(1);
izBlock = round(zBlock./dzBlock);
DP = [-izBlock izBlock];
DP = spdiags(DP, [-1 0],m,m);
DP(1,:) = DP(2,:);
DP = kron(speye(n),DP);
mask = ones(n*m,1);
I = speye(2*m*n);
Ii = I(:,1:m*n);
Id = I(:,1+m*n:end);

% Objective functions
Fid = []; Reg = []; Dual = [];
Fid(1) = 0.5*norm(modelFreq(u,zP,freqVec) - bzf(:))^2;
Reg(1) = muAlpha*TVcalc_isotropic(DP*u(1:m*n),m,n,mask) + ...
    muBeta*TVcalc_isotropic(DP*u(m*n+1:end),m,n,mask);
Dual(1) = 0;
ite  = 0;
error = 1;
tic
while abs(error) > tol && ite < maxIte
    ite = ite + 1;
    
    % Fidelity step
    u = optimNonLinearGN(X,Y, speye(2*m*n),v-w, rho, tol,u);

    % Regularization step
    v = optimAdmmTvDy(Ii,Id,u+w, muAlpha/rho,muBeta/rho ,m,n,tol,mask,DP);

    % Dual update
    w = w + u - v;

    % Loss
    Fid(ite+1) = 1/2*norm( Y - model(u,X) )^2;
    Reg(ite+1) = muAlpha*TVcalc_isotropic(DP*v(1:m*n),m,n,mask) + ...
        muBeta*TVcalc_isotropic(DP*v(m*n+1:end),m,n,mask);
    Dual(ite+1) = norm(u-v);
    % fprintf("Ite: %01i, Fid: %.3f, Reg: %.3f, Dual: %.3f\n",...
    %     ite,Fid(ite+1),Reg(ite+1),Dual(ite+1))
    % fprintf("Ite: %01i, FID+REG = %.3f\n",...
    %     ite,Fid(ite+1)+Reg(ite+1));

    error = Fid(ite+1) + Reg(ite+1) - Fid(ite) - Reg(ite);


    if mod(ite,5)==1
        alphaArr = reshape(DP*u(1:m*n),[m,n]);
        betaArr = reshape(DP*u(1+m*n:end),[m,n]);
        
        estACtv = alphaArr/freq^2;
        estBAtv = 2*(betaArr-1);
        
        
        [Xmesh,Zmesh] = meshgrid(xBlock,zBlock);
        inc = Xmesh.^2 + (Zmesh-centerDepth).^2 < (radiusDisk-1e-3).^2;
        back = Xmesh.^2 + (Zmesh-centerDepth).^2 > (radiusDisk+1e-3).^2;
                
        figure; imagesc(x*1e3,z*1e3,estBAtv); colorbar;
        clim([5 10]);
        title('B/A');
        axis image
        colormap pink; colorbar;
        xlabel('Lateral distance (mm)');
        ylabel('Depth (mm)');
        hold on
        rectangle('Position',[0-radiusDisk,centerDepth-radiusDisk,...
        2*radiusDisk,2*radiusDisk]*1000, 'Curvature',1,...
        'EdgeColor','b', 'LineStyle','--', 'LineWidth',2)
        hold off
        pause(0.1)
        fprintf("Ite: %01i, FID+REG = %.3f\n",...
        ite,Fid(ite+1)+Reg(ite+1));

    end
end
toc


figure('Units','centimeters', 'Position',[5 5 10 15]),
tiledlayout(3,1)
nexttile,
plot(Fid, 'LineWidth',2)
% xlabel('Number of iterations')
ylabel('Fide')
xlim([3 length(Fid)])

nexttile,
plot(Reg, 'LineWidth',2)
% xlabel('Number of iterations')
ylabel('Regularization')
xlim([3 length(Fid)])

nexttile,
plot(Dual, 'LineWidth',2)
xlabel('Number of iterations')
ylabel('Dual step norm')
xlim([3 length(Fid)])

%%
alphaArr = reshape(DP*u(1:m*n),[m,n]);
betaArr = reshape(DP*u(1+m*n:end),[m,n]);
% alphaArr = reshape(u(1:m*n),[m,n]);
% betaArr = reshape(u(1+m*n:end),[m,n]);

estACtv = alphaArr/freq^2;
estBAtv = 2*(betaArr-1);


[Xmesh,Zmesh] = meshgrid(xBlock,zBlock);
inc = Xmesh.^2 + (Zmesh-centerDepth).^2 < (radiusDisk-1e-3).^2;
back = Xmesh.^2 + (Zmesh-centerDepth).^2 > (radiusDisk+1e-3).^2;

fprintf('AC: %.3f +/- %.3f\n', mean(estACtv(:), 'omitnan'), ...
std(estACtv(:), [] ,'omitnan'));
fprintf('B/A inc: %.2f +/- %.2f\n', mean(estBAtv(inc),'omitnan'), ...
std(estBAtv(:), [] ,'omitnan'));
fprintf('B/A back: %.2f +/- %.2f\n', mean(estBAtv(back),'omitnan'), ...
std(estBAtv(:), [] ,'omitnan'));


figure; imagesc(x*1e3,z*1e3,estACtv); colorbar;
clim([0 0.2]);
title('\alpha_0 in \alpha(f) = \alpha_0 \times f^2 dB/cm')
axis image
colormap pink; colorbar;
xlabel('Lateral distance (mm)');
ylabel('Depth (mm)');
hold on
rectangle('Position',[0-radiusDisk,centerDepth-radiusDisk,...
2*radiusDisk,2*radiusDisk]*1000, 'Curvature',1,...
'EdgeColor','b', 'LineStyle','--', 'LineWidth',2)
hold off


figure; imagesc(x*1e3,z*1e3,estBAtv); colorbar;
clim([5 10]);
title('B/A');
axis image
colormap pink; colorbar;
xlabel('Lateral distance (mm)');
ylabel('Depth (mm)');
hold on
rectangle('Position',[0-radiusDisk,centerDepth-radiusDisk,...
2*radiusDisk,2*radiusDisk]*1000, 'Curvature',1,...
'EdgeColor','b', 'LineStyle','--', 'LineWidth',2)
hold off
pause(0.1)


%% Gauss-Newton with TV in step VERSION 2
% Hyperparameters
tol = 1e-3;
muAlpha = 1e-5;
muBeta = 1e-5;
maxIte = 20;

% Auxiliar
dzBlock = zBlock(2)-zBlock(1);
izBlock = round(zBlock./dzBlock);
DP = [-izBlock izBlock];
DP = spdiags(DP, [-1 0],m,m);
DP(1,:) = DP(2,:);
DP = kron(speye(n),DP);

% Initialization
betaIni = 1+(7.5)/2;
alphaIni = 0.1*freq^2; % dB/cm
uC = [alphaIni*ones(n*m,1);betaIni*ones(n*m,1)];
mask = ones(n*m,1);

stepNorm = [];
ite = 1;
tic
while true
    jcb = jacobian(uC,X);
    res = Y - model(uC,X);

    [step,uL] = optimAdmmStepTv(jcb,DP,res,uC,muAlpha,muBeta,m,n,tol,mask);
    uC = uC + step;

    stepNorm(ite) = norm(step);
    fprintf('\tIteration %i,\t step norm: %.3f\n',ite,stepNorm(ite))
    if stepNorm(ite)<tol || ite == maxIte, break; end
    ite = ite + 1;
end
toc
%%

alphaArr = reshape(DP*uC(1:m*n),[m,n]);
betaArr = reshape(DP*uC(1+m*n:end),[m,n]);

estACtv = alphaArr/freq^2;
estBAtv = 2*(betaArr-1);


[Xmesh,Zmesh] = meshgrid(xBlock,zBlock);
inc = Xmesh.^2 + (Zmesh-centerDepth).^2 < (radiusDisk-1e-3).^2;
back = Xmesh.^2 + (Zmesh-centerDepth).^2 > (radiusDisk+1e-3).^2;

fprintf('AC: %.3f +/- %.3f\n', mean(estACtv(:), 'omitnan'), ...
std(estACtv(:), [] ,'omitnan'));
fprintf('B/A inc: %.2f +/- %.2f\n', mean(estBAtv(inc),'omitnan'), ...
std(estBAtv(:), [] ,'omitnan'));
fprintf('B/A back: %.2f +/- %.2f\n', mean(estBAtv(back),'omitnan'), ...
std(estBAtv(:), [] ,'omitnan'));


figure; imagesc(x*1e3,z*1e3,estACtv); colorbar;
clim([0 0.2]);
title('\alpha_0 in \alpha(f) = \alpha_0 \times f^2 dB/cm')
axis image
colormap pink; colorbar;
xlabel('Lateral distance (mm)');
ylabel('Depth (mm)');
hold on
rectangle('Position',[0-radiusDisk,centerDepth-radiusDisk,...
2*radiusDisk,2*radiusDisk]*1000, 'Curvature',1,...
'EdgeColor','b', 'LineStyle','--', 'LineWidth',2)
hold off


figure; imagesc(x*1e3,z*1e3,estBAtv); colorbar;
clim([5 10]);
title('B/A');
axis image
colormap pink; colorbar;
xlabel('Lateral distance (mm)');
ylabel('Depth (mm)');
hold on
rectangle('Position',[0-radiusDisk,centerDepth-radiusDisk,...
2*radiusDisk,2*radiusDisk]*1000, 'Curvature',1,...
'EdgeColor','b', 'LineStyle','--', 'LineWidth',2)
hold off
pause(0.1)
