%% B/A estimator with time domain approach
% by: Andres Coila

clear; close all; clc;
addpath(genpath(pwd))


%% Prepare files for sample and reference
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

mzaxis = movmean(BAeff.mz,4,2);
x = BAeff.lateral';
z = BAeff.axial;
dz = z(2)-z(1);
dx = x(2)-x(1);

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

%% MY VERSION
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
tol = 1e-3;
maxIte = 10;
muAlpha = 100;
muBeta = 1;
beta0 = 1+(10.5)/2;
alpha0 = 0.1*freq^2; % dB/cm -> Np/m

theta = [alpha0*ones(n*m,1);beta0*ones(n*m,1)];
regMatrix = blkdiag(muAlpha*speye(n*m),muBeta*speye(n*m));

stepNorm = [];
ite = 1;
while true
    jcb = jacobian(theta,X);
    res = Y - model(theta,X);
    [step,~] = cgs(jcb'*jcb + regMatrix,jcb'*-res);
    % step = (jcb'*jcb + regMatrix)\(jcb'*-res);
    theta = theta + step;

    stepNorm(ite) = norm(step)/2;
    if stepNorm(ite)<tol || ite == maxIte, break; end
    ite = ite + 1;
end

alphaArr = theta(1:n*m);
betaArr = theta(n*m+1:end);

estAClm = reshape(alphaArr/freq^2,[m,n]);
estBAlm = reshape(2*(betaArr-1),[m,n]);

fprintf('AC: %.3f +/- %.3f\n', mean(estAClm(:), 'omitnan'), ...
    std(estAClm(:), [] ,'omitnan'));
fprintf('B/A: %.2f +/- %.2f\n', mean(estBAlm(:), 'omitnan'), ...
    std(estBAlm(:), [] ,'omitnan'));

figure('Units','centimeters', 'Position',[5 5 10 5]), 
plot(stepNorm, 'LineWidth',2)
xlabel('Number of iterations')
ylabel('Step norm')
axis tight

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
muLocal = 0.03;
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

%% Gauss-Newton with TV in step VERSION 2
% Hyperparameters
tol = 1e-3;
muAlpha = 1e-5;
muBeta = 1e-5;

% Auxiliar
dzBlock = zBlock(2)-zBlock(1);
izBlock = round(zBlock./dzBlock);
DP = [-izBlock izBlock];
DP = spdiags(DP, [-1 0],m,m);
DP(1,:) = DP(2,:);
DP = kron(speye(n),DP);
mask = ones(n*m,1);

% Initialization
beta0 = 1+(7.5)/2;
alpha0 = 0.08*freq^2; % dB/cm
uC = [alpha0*ones(n*m,1);beta0*ones(n*m,1)];
jcb = jacobian(uC,X);
res = Y - model(uC,X);

% v = [alpha0*ones(n*m,1);beta0*ones(n*m,1)];
% w = zeros(size(u));
% alphaCk = theta(1:m*n);
% betaCk = theta(m*n+1:end);
[step,uL] = optimAdmmStepTv(jcb,DP,res,uC,muAlpha,muBeta,m,n,tol,mask);

%% Many iterations
beta0 = 1+(7.5)/2;
alpha0 = 0.1*freq^2; % dB/cm
uC = [alpha0*ones(n*m,1);beta0*ones(n*m,1)];
mask = ones(n*m,1);
tol = 1e-3;
maxIte = 20;
muAlpha = 1e-2;
muBeta = 1e-2;

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
figure('Units','centimeters', 'Position',[5 5 10 5]), 
plot(stepNorm, 'LineWidth',2)
xlabel('Number of iterations')
ylabel('Step norm')
axis tight

alphaArr = DP*uC(1:m*n);
betaArr = DP*uC(1+m*n:end);
estACtv = reshape(alphaArr/freq^2,[m,n]);
estBAtv = reshape(2*(betaArr-1),[m,n]);

fprintf('AC: %.3f +/- %.3f\n', mean(estACtv(:), 'omitnan'), ...
std(estACtv(:), [] ,'omitnan'));
fprintf('B/A: %.2f +/- %.2f\n', mean(estBAtv(:), 'omitnan'), ...
std(estBAtv(:), [] ,'omitnan'));

figure; imagesc(x*1e3,z*1e3,estACtv); colorbar;
clim([0 0.2]);
title('\alpha_0 in \alpha(f) = \alpha_0 \times f^2 dB/cm')
font = 20;
axis image
colormap turbo; colorbar;
set(gca,'fontsize',font)
xlabel('Lateral distance (mm)');
ylabel('Depth (mm)');
set(gca,'FontSize',font);
hold on
rectangle('Position',[0-radiusDisk,centerDepth-radiusDisk,...
2*radiusDisk,2*radiusDisk]*1000, 'Curvature',1,...
'EdgeColor','b', 'LineStyle','--', 'LineWidth',2)
hold off


figure; imagesc(x*1e3,z*1e3,estBAtv); colorbar;
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

%% Gauss-Newton with TV in step
% Hyperparameters
tol = 1e-3;
muAlpha = 1e1;
muBeta = 1e1;

% Auxiliar
dzBlock = zBlock(2)-zBlock(1);
izBlock = round(zBlock./dzBlock);
DP = [-izBlock izBlock];
DP = spdiags(DP, [-1 0],m,m);
DP(1,:) = DP(2,:);
DP = kron(speye(n),DP);
mask = ones(n*m,1);
jcb = jacobian(theta,X);
jcbAlpha = jcb(:,1:n*m);
jcbBeta = jcb(:,n*m+1:end);
res = Y - model(theta,X);

% Initialization
beta0 = 1+(10.5)/2;
alpha0 = 0.08*freq^2; % dB/cm -> Np/m
u = [alpha0*ones(n*m,1);beta0*ones(n*m,1)];
v = [alpha0*ones(n*m,1);beta0*ones(n*m,1)];
w = zeros(size(u));
alphaCk = theta(1:m*n);
betaCk = theta(m*n+1:end);
[A,B,C,D] = optimStepTv(jcbAlpha,jcbBeta,DP,res,alphaCk,betaCk,...
    muAlpha,muBeta,m,n,tol,mask);  

%%
beta0 = 1+(10.5)/2;
alpha0 = 0.12*5^2/8.686*100; % dB/cm -> Np/m

theta = [alpha0*ones(n*m,1);beta0*ones(n*m,1)];
mask = ones(n*m,1);
tol = 1e-3;
maxIte = 10;
muAlpha = 0.1; % 10 funciona bien
muBeta = 0.1;

stepNorm = [];
ite = 1;
tic
while true
    jcb = jacobian(theta,X);
    jcbAlpha = jcb(:,1:n*m);
    jcbBeta = jcb(:,n*m+1:end);
    res = Y - model(theta,X);

    alphaCk = theta(1:m*n);
    betaCk = theta(m*n+1:end);
    [A,B,C,D] = optimStepTv(jcbAlpha,jcbBeta,DP,res,alphaCk,betaCk,...
        muAlpha,muBeta,m,n,tol,mask);  

    theta = theta + [A;B];

    stepNorm(ite) = norm([A;B])/2;
    fprintf('\tIteration %i,\t step norm: %.3f\n',ite,stepNorm(ite))
    if stepNorm(ite)<tol || ite == maxIte, break; end
    ite = ite + 1;


end
toc
%%
figure('Units','centimeters', 'Position',[5 5 10 5]), 
plot(stepNorm, 'LineWidth',2)
xlabel('Number of iterations')
ylabel('Step norm')
axis tight


alphaArr = C;
betaArr = D;
alpha_dB = alphaArr/5^2/100*8.686; % Np/cm/MH^2 : 0.1

estACtv = reshape(alpha_dB,[m,n]);
estBAtv = reshape(2*(betaArr-1),[m,n]);

fprintf('AC: %.3f +/- %.3f\n', mean(estACtv(:), 'omitnan'), ...
std(estACtv(:), [] ,'omitnan'));
fprintf('B/A: %.2f +/- %.2f\n', mean(estBAtv(:), 'omitnan'), ...
std(estBAtv(:), [] ,'omitnan'));

figure; imagesc(x*1e3,z*1e3,estACtv); colorbar;
clim([0 0.2]);
title('\alpha_0 in \alpha(f) = \alpha_0 \times f^2 dB/cm')
font = 20;
axis image
colormap jet; colorbar;
set(gca,'fontsize',font)
xlabel('Lateral distance (mm)');
ylabel('Depth (mm)');
set(gca,'FontSize',font);
hold on
rectangle('Position',[0-radiusDisk,centerDepth-radiusDisk,...
2*radiusDisk,2*radiusDisk]*1000, 'Curvature',1,...
'EdgeColor','b', 'LineStyle','--', 'LineWidth',2)
hold off


figure; imagesc(x*1e3,z*1e3,estBAtv); colorbar;
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

%% Other approach: variable change (DOES NOT WORK)
% % Hyperparameters
% tol = 1e-3;
% muAlpha = 1e-2;
% muBeta = 1e-4;
% 
% % Auxiliar
% dzBlock = zBlock(2)-zBlock(1);
% izBlock = round(zBlock./dzBlock);
% DP = [-izBlock izBlock];
% DP = spdiags(DP, [-1 0],m,m);
% DP(1,1) = DP(1,1);
% DP = kron(speye(n),DP);
% P = inv(DP);
% mask = ones(length(Y),1);
% jcb = jacobian(theta,X);
% jcbAlpha = jcb(:,1:n*m);
% jcbBeta = jcb(:,n*m+1:end);
% res = Y - model(theta,X);
% 
% % Initialization
% beta0 = 1+(10.5)/2;
% alpha0 = 0.08*freq^2; % dB/cm -> Np/m
% u = [alpha0*ones(n*m,1);beta0*ones(n*m,1)];
% v = [alpha0*ones(n*m,1);beta0*ones(n*m,1)];
% w = zeros(size(u));
% alphaCk = theta(1:m*n);
% betaCk = theta(m*n+1:end);
% 
% [localAlpha,localBeta] = AlterOpti_ADMM(jcbAlpha*P,jcbBeta*P,...
%     -res+jcbAlpha*alphaCk+jcbBeta*betaCk,...
%     muAlpha,muBeta,m,n,tol,mask);

