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
beta0 = 1+(10.5)/2;
alpha0 = 0.10*5^2/8.686*100; % dB/cm -> Np/m
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
maxIte = 50;
muAlpha = 1;
muBeta = 1;

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
alpha_dB = alphaArr/5^2/100*8.686; % Np/cm/MH^2 : 0.1

estAClm = reshape(alpha_dB,[m,n]);
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

%% Gauss-Newton with TV in step
% Iterating
beta0 = 1+(10.5)/2;
alpha0 = 0.1*5^2/8.686*100; % dB/cm -> Np/m

theta = [alpha0*ones(n*m,1);beta0*ones(n*m,1)];
mask = ones(n*m*np,1);
tol = 1e-3;
maxIte = 20;
muAlpha = 1;
muBeta = 1;

stepNorm = [];
ite = 1;
tic
while true
    jcb = jacobian(theta,X);
    jcbAlpha = jcb(:,1:n*m);
    jcbBeta = jcb(:,n*m+1:end);
    res = Y - model(theta,X);
    [stepAlpha,stepBeta] = AlterOpti_ADMM(jcbAlpha,jcbBeta,-res,...
        muAlpha,muBeta,m,n,tol,mask);
    theta = theta + [stepAlpha;stepBeta];

    stepNorm(ite) = norm([stepAlpha;stepBeta])/2;
    fprintf('\tIteration %i,\t step norm: %.3f\n',ite,stepNorm(ite))
    if stepNorm(ite)<tol || ite == maxIte, break; end
    ite = ite + 1;


end
toc

figure('Units','centimeters', 'Position',[5 5 10 5]), 
plot(stepNorm, 'LineWidth',2)
xlabel('Number of iterations')
ylabel('Step norm')
axis tight


alphaArr = theta(1:n*m);
betaArr = theta(n*m+1:end);
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
colormap pink; colorbar;
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

%% Gauss-Newton with regularization objective
dzBlock = zBlock(2)-zBlock(1);
izBlock = round(zBlock./dzBlock);
DP = [-izBlock izBlock];
DP = spdiags(DP, [-1 0],m,m);
% DP(1,:) = DP(2,:);
DP(1,1) = 1;
% figure, imagesc(DP)
% axis image
DP = kron(speye(n),DP);

tol = 1e-3;
maxIte = 1;
muAlpha = 1;
muBeta = 1;

theta = [alpha0*ones(n*m,1);beta0*ones(n*m,1)];

stepNorm = [];
ite = 1;
while true
    jcb = jacobian(theta,X);
    res = Y - model(theta,X);

    alphaK = theta(1:n*m);
    betaK = theta(n*m+1:end);
    regAlpha = muAlpha*DP*gradTV(DP*alphaK,m,n);
    regBeta = muBeta*DP*gradTV(DP*betaK,m,n);
    % regAlpha = muAlpha*gradTV(alphaK,m,n);
    % regBeta = muBeta*gradTV(betaK,m,n);

    [stepApprox,~] = pcg(jcb'*jcb,-jcb'*res);
    [step,~] = pcg(jcb'*jcb,-jcb'*res - [regAlpha;regBeta],tol,200,[],[],stepApprox);
     % [step,~] = pcg(jcb'*jcb,-jcb'*res - [regAlpha;regBeta]);
    % step = (jcb'*jcb)\(jcb'*-res - [regAlpha;regBeta]);
    theta = theta + step;

    stepNorm(ite) = norm(step)/2;
    fprintf('\tIteration %i,\t step norm: %.3f\n',ite,stepNorm(ite))
    if stepNorm(ite)<tol || ite == maxIte, break; end
    ite = ite + 1;
end

alphaArr = theta(1:n*m);
betaArr = theta(n*m+1:end);
alpha_dB = alphaArr/5^2/100*8.686; % Np/cm/MH^2 : 0.1

estAClm = reshape(alpha_dB,[m,n]);
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

%% Gauss-Newton with TV(DP*u) in step
beta0 = 1+(10.5)/2;
alpha0 = 0.10*5^2/8.686*100; % dB/cm -> Np/m

dzBlock = zBlock(2)-zBlock(1);
izBlock = round(zBlock./dzBlock);
DP = [-izBlock izBlock];
DP = spdiags(DP, [-1 0],m,m);
% DP(1,:) = DP(2,:);
DP(1,1) = 1;
% figure, imagesc(DP)
% axis image
DP = kron(speye(n),DP);


% Iterating
theta = [alpha0*ones(n*m,1);beta0*ones(n*m,1)];
mask = ones(n*m,1);
tol = 1e-3;
maxIte = 5;
muAlpha = 10;
muBeta = 0.1;

stepNorm = [];
ite = 1;
tic
while true
    jcb = jacobian(theta,X);
    jcbAlpha = jcb(:,1:n*m);
    jcbBeta = jcb(:,n*m+1:end);

    res = Y - model(theta,X);
    [stepAlpha,stepBeta] = optimAdmmTvDy(jcbAlpha,jcbBeta,-res,...
        muAlpha,muBeta,m,n,tol,mask,DP);
    theta = theta + [stepAlpha;stepBeta];

    stepNorm(ite) = norm([stepAlpha;stepBeta])/2;
    fprintf('\tIteration %i,\t step norm: %.3f\n',ite,stepNorm(ite))
    if stepNorm(ite)<tol || ite == maxIte, break; end
    ite = ite + 1;


end
toc

figure('Units','centimeters', 'Position',[5 5 10 5]), 
plot(stepNorm, 'LineWidth',2)
xlabel('Number of iterations')
ylabel('Step norm')
axis tight


alphaArr = theta(1:n*m);
betaArr = theta(n*m+1:end);
alpha_dB = alphaArr/5^2/100*8.686; % Np/cm/MH^2 : 0.1

estACtv = reshape(alpha_dB,[m,n]);
estBAtv = reshape(2*(betaArr-1),[m,n]);

fprintf('AC: %.3f +/- %.3f\n', mean(estACtv(:), 'omitnan'), ...
std(estACtv(:), [] ,'omitnan'));
fprintf('B/A: %.2f +/- %.2f\n', mean(estBAtv(:), 'omitnan'), ...
std(estBAtv(:), [] ,'omitnan'));

figure; imagesc(x*1e3,z*1e3,estACtv); colorbar;
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
hold on
rectangle('Position',[0-radiusDisk,centerDepth-radiusDisk,...
2*radiusDisk,2*radiusDisk]*1000, 'Curvature',1,...
'EdgeColor','b', 'LineStyle','--', 'LineWidth',2)
hold off
pause(0.5)


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

%% Inversion
dzBlock = zBlock(2)-zBlock(1);
izBlock = round(zBlock./dzBlock);
DP = [-izBlock izBlock];
DP = spdiags(DP, [-1 0],m,m);
% DP(1,:) = DP(2,:);
DP(1,1) = 1;
DP = kron(speye(n),DP);

estAClocal = reshape(DP*estACtv(:),[m,n]);
estBAlocal = reshape(DP*estBAtv(:),[m,n]);

figure; imagesc(x*1e3,z*1e3,estAClocal); colorbar;
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
hold on
rectangle('Position',[0-radiusDisk,centerDepth-radiusDisk,...
2*radiusDisk,2*radiusDisk]*1000, 'Curvature',1,...
'EdgeColor','b', 'LineStyle','--', 'LineWidth',2)
hold off
pause(0.5)


figure; imagesc(x*1e3,z*1e3,estBAlocal); colorbar;
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

%% GN with TV(DP*u) in step and change of variable
dzBlock = zBlock(2)-zBlock(1);
izBlock = round(zBlock./dzBlock);
DP = [-izBlock izBlock];
DP = spdiags(DP, [-1 0],m,m);
% DP(1,:) = DP(2,:);
DP(1,1) = 1;
figure, imagesc(DP)
axis image
DP = kron(speye(n),DP);
P = inv(DP);

% Iterating
theta = [alpha0*ones(n*m,1);beta0*ones(n*m,1)];
mask = ones(n*m*np,1);
tol = 1e-3;
maxIte = 2;
muAlpha = 10;
muBeta = 1;

stepNorm = [];
ite = 1;
while true
    jcb = jacobian(theta,X);
    jcbAlpha = jcb(:,1:n*m);
    jcbBeta = jcb(:,n*m+1:end);

    alphaCum = theta(1:n*m);
    betaCum = theta(n*m+1:end);

    res = Y - model(theta,X);
    tic
    [localAlpha,localBeta] = AlterOpti_ADMM(jcbAlpha*P,jcbBeta*P,...
        -res+jcbAlpha*alphaCum+jcbBeta*betaCum,...
        muAlpha,muBeta,m,n,tol,mask);
    toc
    nextAlpha = P*localAlpha;
    nextBeta = P*localBeta;
    theta = [nextAlpha;nextBeta];

    loss(ite) = norm(res).^2/2 + TVcalc_isotropic(localAlpha,m,n,ones(m,n)) + ...
        + TVcalc_isotropic(localBeta,m,n,ones(m,n));
    % fprintf('\tIteration %i,\t step norm: %.3f\n',ite,stepNorm(ite))
    if ite == maxIte, break; end
    ite = ite + 1;

    % Plot
    alphaArr = localAlpha;
    betaArr = localBeta;
    alpha_dB = alphaArr/5^2/100*8.686; % Np/cm/MH^2 : 0.1
    estACtv = reshape(alpha_dB,[m,n]);
    estBAtv = reshape(2*(betaArr-1),[m,n]);

    figure; imagesc(x*1e3,z*1e3,estACtv); colorbar;
    clim([0 0.2]);
    title('\alpha_0 in \alpha(f) = \alpha_0 \times f^2 dB/cm')
    font = 20;
    axis image
    colormap pink; colorbar;
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

end


figure('Units','centimeters', 'Position',[5 5 10 5]), 
plot(stepNorm, 'LineWidth',2)
xlabel('Number of iterations')
ylabel('Step norm')
axis tight


% alphaArr = theta(1:n*m);
% betaArr = theta(n*m+1:end);


%% Utility functions

function jcb = jacobian(theta,x)
mn = size(theta,1)/2; % number of points in image
p = length(x)/mn;

alphaArr = theta(1:mn)';
betaArr = theta(mn+1:end)';
zBA = reshape(x,[p,mn]);

dMLfda = -(betaArr).*((2*zBA.*alphaArr+1).*...
    exp(-2*alphaArr.*zBA) - 1)./(alphaArr.^2.*zBA);
dMLfdb = -(1./zBA./alphaArr).*(1 - exp(-2*alphaArr.*zBA));

% Indices for sparse matrix
I = reshape(1:mn*p,[mn,p]); % Row indices
J = (1:mn).*ones(p,1);

jcb = [sparse(I,J,dMLfda) sparse(I,J,dMLfdb)];
end

function mdl = model(theta,x)
mn = size(theta,1)/2; % number of points in image
p = length(x)/mn;

alphaArr = theta(1:mn)';
betaArr = theta(mn+1:end)';
zBA = reshape(x,[p,mn]);

mdl = (betaArr./zBA./alphaArr).*(1 - exp(-2*alphaArr.*zBA)) ;
mdl = mdl(:);
end

function gTV = gradTV(im,m,n)
    im = reshape(im,m,n);
    [GX,GY] = gradient(im);
    eps = 0.3;
    Gmod = sqrt(GX.^2 + GY.^2 + eps);
    gTV = -divergence(GX./Gmod,GY./Gmod);
    gTV = gTV(:);
end