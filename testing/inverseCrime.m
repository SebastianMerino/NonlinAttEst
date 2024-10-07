%% B/A estimator with time domain approach
% by: Andres Coila

clear; close all; clc;
addpath(genpath(pwd))


% Prepare files for sample and reference
baseDir = ['C:\Users\sebas\Documents\MATLAB\totalvarsimul_AC_BA\' ...
    'BA_AC_joint\maps'];
% fileSam = 'rf_fnum3_PWNE_samBA6_att0p10f2_nc10_400kPa';
% fileRef = 'rf_fnum3_PWNE_refBA6_att10f2_nc10_400kPa';
fileSam = 'rf_fnum3_SCOMP5MHz_nc10_0p10f2_saminc400_doubleangle720';
% fileSam = 'rf_fnum3_SCOMP5MHz_nc10_0p10f2_sam400_doubleangle720';
fileRef = 'rf_fnum3_SCOMP5MHz_nc10_0p10f2_ref400_doubleangle720';

filename = ['FULLMAPv60_',fileSam,'_',fileRef,'.mat'];
load(fullfile(baseDir,filename));

%% Getting system
freq = 5;
zlim = [5 55]*1e-3;
nwl = 20;
radiusDisk = (9)*1e-3;
centerDepth = 22.5e-3;

x = BAeff.lateral';
z = BAeff.axial;
dz = z(2)-z(1);
dx = x(2)-x(1);

%% Generating local maps
NptodB = db(exp(1));
[Xmesh,Zmesh] = meshgrid(x,z);
inc = Xmesh.^2 + (Zmesh-centerDepth).^2 <= (radiusDisk).^2;
betaL = ones(size(Xmesh))*(1+6/2);
% betaL(inc) = (1+9/2);
alphaL = ones(size(Xmesh))*0.1*freq^2;

figure('Units','centimeters', 'Position',[5 5 20 10])
font = 9;
tiledlayout(1,2)
t1 = nexttile;
imagesc(x,z,alphaL/freq^2); colorbar;
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
imagesc(x*1e2,z*1e2,alphaC/freq^2); colorbar;
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

alphaCnp = alphaC/NptodB*100; % dB/cm -> Np/m
mzaxis = betaC.*(1 - exp(-2*alphaCnp.*Zmesh) )./alphaCnp./Zmesh + ...
    0.0*randn(size(Xmesh));
%% Getting system
wl = 1540/1e6/freq;
blocksize = nwl*wl; % Change to test different sizes
nz = floor(blocksize/dz);
nx = floor(blocksize/dx);

zori = z;
mzaxisori = mzaxis;
[~,idzmin] = min(abs(zlim(1)-z));
[~,idzmax] = min(abs(zlim(2)-z));
z = zori(idzmin:idzmax);
mzaxis = mzaxisori(idzmin:idzmax,:);
L1   = size(mzaxis,1);
L2   = size(mzaxis,2);

param.overlap_pct = 0.1;
param.nz = nz;
param.overlap = round((1-param.overlap_pct)*param.nz);
param.nooverlap = nz - param.overlap;
zini(1)  = 1;
dz = 1e-3;

while zini(end) < L1 + 1 - nz - param.nooverlap
    zini(end+1) = zini(end)+param.nooverlap;
end

zfin = zini + nz - 1;
f0 = freq;
m = length(zfin);
n = L2;

% Constructing arrays
zblock = (z(zini(1)):1e-3:z(zfin(1))); % m
np = length(zblock); % Number of points per block
mzBA = zeros(np,m*n);
zBA = zeros(np,m*n); % zBA
for pixj = 1:n
    for pixi = 1:m
        zblock = (z(zini(pixi)):1e-3:z(zfin(pixi))); % m
        for zi = 1:length(zblock)
            [~, idz] = min(abs(zblock(zi)-z));
            zBA(zi,pixi+(pixj-1)*m) = z(idz);
            mzBA(zi,pixi+(pixj-1)*m) = mzaxis(idz,pixj);
            % X((pixj-1)*np+(pixi-1)*n*np+zi) = z(idz);
            % Y((pixj-1)*np+(pixi-1)*n*np+zi) = mzaxis(idz,pixj);
        end
    end
end
zBlock = (reshape(zBA(round(np/2),:),m,n));
zBlock = zBlock(:,1);
xBlock = x;
X = zBA(:); % zBA
Y = mzBA(:); % mzBA

%% Gauss-Newton with LM
tol = 1e-5;
maxIte = 400;
muAlpha = 0;
muBeta = 0;
% muAlpha = 0; muBeta = 0;
beta0 = 1+(10.5)/2;
alpha0 = 0.15*freq^2; % dB/cm -> Np/m

theta = [alpha0*ones(n*m,1);beta0*ones(n*m,1)];
regMatrix = blkdiag(muAlpha*speye(n*m),muBeta*speye(n*m));

loss = [];
loss(1) = 0.5*norm(Y - model(theta,X))^2;

ite = 1;
tic 
while true
    jcb = jacobian(theta,X);
    res = (Y - model(theta,X));
    [step,~] = pcg(jcb'*jcb + regMatrix,jcb'*-res, tol, 200);
    % step = (jcb'*jcb + regMatrix)\(jcb'*-res);
    theta = theta + step;

    loss(ite+1) = 0.5*norm(Y - model(theta,X))^2;
    if abs(loss(ite+1) - loss(ite))<tol || ite == maxIte, break; end
    ite = ite + 1;
end
toc
%%
alphaArr = theta(1:n*m);
betaArr = theta(n*m+1:end);

estAClm = reshape(alphaArr/freq^2,[m,n]);
estBAlm = reshape(2*(betaArr-1),[m,n]);

fprintf('AC: %.3f +/- %.3f\n', mean(estAClm(:), 'omitnan'), ...
    std(estAClm(:), [] ,'omitnan'));
fprintf('B/A: %.2f +/- %.2f\n', mean(estBAlm(:), 'omitnan'), ...
    std(estBAlm(:), [] ,'omitnan'));

figure('Units','centimeters', 'Position',[5 5 10 5]), 
plot(loss, 'LineWidth',2)
xlim([1 length(loss)])
xlabel('Number of iterations')
ylabel('Loss')

figure; imagesc(x*1e3,z*1e3,estAClm); colorbar;
clim([0 0.2]);
%title('ACS (dB/cm/MHz)');
title('\alpha_0 in \alpha(f) = \alpha_0 \times f^2 dB/cm')
font = 20;
axis image
colormap pink; colorbar;
%clim([0 12])
%font = 20;
set(gca,'fontsize',font)
xlabel('Lateral distance (mm)');
ylabel('Depth (mm)');
%title("\bf BoA map");
set(gca,'FontSize',font);
%ylim([5 55])
hold on
rectangle('Position',[0-radiusDisk,centerDepth-radiusDisk,...
    2*radiusDisk,2*radiusDisk]*1000, 'Curvature',1,...
    'EdgeColor','b', 'LineStyle','--', 'LineWidth',2)
hold off


figure; imagesc(x*1e3,z*1e3,estBAlm); colorbar;
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
hold on
rectangle('Position',[0-radiusDisk,centerDepth-radiusDisk,...
    2*radiusDisk,2*radiusDisk]*1000, 'Curvature',1,...
    'EdgeColor','b', 'LineStyle','--', 'LineWidth',2)
hold off
pause(0.5)

%% Inversion, Coila's implementation
muLocal = 0.1;
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
beta0 = 1+(7.5)/2;
alpha0 = 0.1*freq^2; % dB/cm
u = [alpha0*ones(n*m,1);beta0*ones(n*m,1)];
v = [alpha0*ones(n*m,1);beta0*ones(n*m,1)];
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
Fid(1) = 1/2*norm( Y - model(u,X) )^2;
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
    % v = optimAdmmMixed(Ii,Id,u+w, muAlpha/rho,muBeta/rho ,m,n,tol,mask,DP);
    % [vi,vd] = AlterOpti_ADMM(Ii,Id,u+w, muAlpha/rho,muBeta/rho ,m,n,tol,ones(2*m*n,1));
    % v = [vi;vd];

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
        
        % fprintf('AC: %.3f +/- %.3f\n', mean(estACtv(:), 'omitnan'), ...
        % std(estACtv(:), [] ,'omitnan'));
        % fprintf('B/A inc: %.2f +/- %.2f\n', mean(estBAtv(inc),'omitnan'), ...
        % std(estBAtv(:), [] ,'omitnan'));
        % fprintf('B/A back: %.2f +/- %.2f\n', mean(estBAtv(back),'omitnan'), ...
        % std(estBAtv(:), [] ,'omitnan'));
        
        % figure; imagesc(x*1e3,z*1e3,estACtv); colorbar;
        % clim([0 0.2]);
        % title('\alpha_0 in \alpha(f) = \alpha_0 \times f^2 dB/cm')
        % axis image
        % colormap pink; colorbar;
        % xlabel('Lateral distance (mm)');
        % ylabel('Depth (mm)');
        % hold on
        % rectangle('Position',[0-radiusDisk,centerDepth-radiusDisk,...
        % 2*radiusDisk,2*radiusDisk]*1000, 'Curvature',1,...
        % 'EdgeColor','b', 'LineStyle','--', 'LineWidth',2)
        % hold off
        
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
beta0 = 1+(7.5)/2;
alpha0 = 0.1*freq^2; % dB/cm
uC = [alpha0*ones(n*m,1);beta0*ones(n*m,1)];
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
