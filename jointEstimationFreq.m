% My implementation of the new method with frequency compounding
% Used for testing several frequencies

clear; close all; clc;
addpath(genpath(pwd))

baseDir = 'C:\Users\sebas\Documents\Data\Nonlinearity\attInc';
% baseDir = 'C:\Users\smerino.C084288\Documents\MATLAB\BA_AC_joint\rfdata';

% fileSam = 'rf_fnum3_SCOMP5MHz_nc10_0p10f2_saminc400_doubleangle720';
% fileRef = 'rf_fnum3_SCOMP5MHz_nc10_0p10f2_ref400_doubleangle720';

% fileSam = 'rf_baBack6_baInc9_att0p1.mat';
% fileSam = 'rf_ba9_attBack0p1_attInc0p18.mat';
% fileSam = 'rf_baBack6_baInc9_attBack0p1_attInc0p18.mat';
% fileRef = 'rf_ba8_att0p12_ref.mat';

% fileSam = 'RFfn3_PWNE_samBA9_att0p10f2_nc10_400kPa';
fileSam = 'RFfn3_PWNE5MHz_samincBA6inc12_att0p1f2inc0p14_nc10_400kPa';
% fileSam = 'RFfn3_PWNE5MHz_samincBA6inc12_att0p1f2inc0p10_nc10_400kPa';
% fileRef = 'RFfn3_PWNE_refBA6_att8f2_nc10_400kPa';
% fileRef = 'RFfn3_PWNE_refBA6_att10f2_nc10_400kPa';
fileRef = 'RFfn3_PWNE_refBA6_att12f2_nc10_400kPa';

% Auxiliar variables
NptodB = 20*log10(exp(1));
% radiusDisk = (9)*1e-3;
% centerDepth = 22.5e-3;
radiusDisk = (9)*1e-3;
centerDepth = 22.5e-3;

baRange = [5 13];
attRange = [0.05,0.2];

%% Preparing data
% Known variables
medium.v = 5;
medium.betaR = 1 + 6/2;             
medium.alphaR = 0.12/NptodB*100; % alpha0 in dB/100/MHz2

% Sample
sample = load(fullfile(baseDir,fileSam));
medium.z = sample.z';
medium.x = sample.x;
medium.fs = sample.fs;
medium.rfL = sample.rf1(:,:,1);
medium.rfH = sample.rf2(:,:,1);
clear sample

% Reference
ref = load(fullfile(baseDir,fileRef));
medium.rfLR = ref.rf1(:,:,1);
medium.rfHR = ref.rf2(:,:,1);
clear ref

% Filtering parameters
filterParams.freqC = 5;
filterParams.freqTol = 2;
filterParams.nCycles = 10; % Number of cycles of the initial filter

% Subsampling parameters
wl = 1540/5/1e6; % Mean central frequency
blockParams.blockSize = [15 15]*wl; 
blockParams.overlap = 0.8;
blockParams.zlim = [0.5; 5]/100;
blockParams.xlim = [-2.5; 2.5]/100;

freqVec = [4,5,6]; % FRECUENCIES FOR FILTERING

%% Bmode
bmode = db(hilbert(medium.rfL(200:end-600,:,1)));
bmode = bmode - max(bmode (:));
figure,
imagesc(medium.x*100,medium.z(200:end-600)*100,bmode, [-40 0])
colormap gray
axis image
colorbar
title('B-mode')
%% Getting B/A
bzf = [];
for iFreq = 1:length(freqVec)
    filterParams.freqC = freqVec(iFreq);
    [bz,xP,zP] = getMeasurements(medium,filterParams,blockParams);
    bzf(:,:,iFreq) = bz;
end
%% Measurements
figure('Units','centimeters', 'Position',[5 5 20 6]),
tiledlayout(1,3)
nexttile,
imagesc(xP*100,zP*100,bzf(:,:,1), [2 7])
title("Measurements b(z,f)")
xlabel('Lateral [cm]')
ylabel('Depth [cm]')
axis image
colorbar
colormap parula

nexttile,
imagesc(xP*100,zP*100,bzf(:,:,2), [2 7])
title("Measurements b(z,f)")
xlabel('Lateral [cm]')
ylabel('Depth [cm]')
axis image
colorbar
colormap parula

nexttile,
imagesc(xP*100,zP*100,bzf(:,:,3), [2 7])
title("Measurements b(z,f)")
xlabel('Lateral [cm]')
ylabel('Depth [cm]')
axis image
colorbar
colormap parula

%% Initialization
[m,n,p] = size(bzf);

[Xmesh,Zmesh] = meshgrid(xP,zP);
inc = Xmesh.^2 + (Zmesh-centerDepth).^2 <= (radiusDisk).^2;
alphaL = ones(size(Xmesh))*0.1;
alphaL(inc) = 0.14;
betaL = ones(size(Xmesh))*(1+7.5/2);

izBlock = 1:length(zP);
P = sparse(tril(ones(m)));
P = P./izBlock';

betaC = P*(betaL);
alphaC =  P*(alphaL);

u0 = [alphaC(:)*100/NptodB;betaC(:)];

figure,
tiledlayout(1,2)
nexttile,
imagesc(xP*1e3,zP*1e3,alphaC)
clim(attRange);
title('Init \alpha_0 in \alpha(f) = \alpha_0 \times f^2')
axis image
colormap turbo; colorbar;
xlabel('Lateral distance (mm)');
ylabel('Depth (mm)');

t2 = nexttile;
imagesc(xP*1e3,zP*1e3,(betaC-1)*2)
clim(baRange);
title('Init B/A');
axis image
colormap(t2,pink); colorbar;
xlabel('Lateral distance (mm)');
ylabel('Depth (mm)');
pause(0.5)

%% Gauss-Newton with LM
[m,n,p] = size(bzf);
tol = 1e-3;
maxIte = 400;
muAlpha = 0;
muBeta = 0;
regMatrix = blkdiag(muAlpha*speye(n*m),muBeta*speye(n*m));

u = u0;
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

    alphaArr = u(1:n*m);
    betaArr = u(n*m+1:end);

    estAClm = reshape(alphaArr*NptodB/100,[m,n]);
    estBAlm = reshape(2*(betaArr-1),[m,n]);
    
    if ite<5
        figure,
        tiledlayout(1,2)
        nexttile,
        imagesc(xP*1e3,zP*1e3,estAClm)
        clim(attRange);
        title('\alpha_0 in \alpha(f) = \alpha_0 \times f^2')
        axis image
        colormap turbo; colorbar;
        xlabel('Lateral distance (mm)');
        ylabel('Depth (mm)');
        
        t2 = nexttile;
        imagesc(xP*1e3,zP*1e3,estBAlm)
        clim(baRange);
        title('B/A');
        axis image
        colormap(t2,pink); colorbar;
        xlabel('Lateral distance (mm)');
        ylabel('Depth (mm)');
        pause(0.5)
    end

end
toc

alphaArr = u(1:n*m);
betaArr = u(n*m+1:end);

estAClm = reshape(alphaArr*NptodB/100,[m,n]);
estBAlm = reshape(2*(betaArr-1),[m,n]);

fprintf('AC: %.3f +/- %.3f\n', mean(estAClm(:), 'omitnan'), ...
std(estAClm(:), [] ,'omitnan'));
fprintf('B/A : %.2f +/- %.2f\n', mean(estBAlm(:),'omitnan'), ...
std(estBAlm(:), [] ,'omitnan'));

figure('Units','centimeters', 'Position',[5 5 10 5]), 
plot(loss, 'LineWidth',2)
xlim([2 length(loss)])
xlabel('Number of iterations')
ylabel('Loss')


figure; imagesc(xP*1e3,zP*1e3,estAClm); colorbar;
clim(attRange);
title('\alpha_0 in \alpha(f) = \alpha_0 \times f^2 dB/100')
axis image
colormap turbo; colorbar;
xlabel('Lateral distance (mm)');
ylabel('Depth (mm)');

figure; imagesc(xP*1e3,zP*1e3,estBAlm); colorbar;
clim(baRange);
title('B/A');
axis image
colormap pink; colorbar;
xlabel('Lateral distance (mm)');
ylabel('Depth (mm)');
pause(0.5)

%% Local maps with regularization
muLocal = 0.1;
dzP = zP(2)-zP(1);
izP = round(zP./dzP);

factorq = izP(1)./izP ;
estBAcum = estBAlm - estBAlm(1,:).*factorq; 
% estBAcum = estBAtv - estBAtv(1,:).*factorq;
estBAcum = estBAcum(2:end,:);


P = sparse(tril(ones(m-1)));
P = P./izP(2:end);
P = kron(speye(n),P);
estBAinst = IRLS_TV(estBAcum(:),P,muLocal,m-1,n,tol,[],ones((m-1)*n,1));
estBAinst = reshape(estBAinst,m-1,n);

figure; imagesc(xP*1e3,zP(2:end)*1e3,estBAinst); colorbar;
clim(baRange);
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

[Xmesh,Zmesh] = meshgrid(xP,zP(2:end));
inc = Xmesh.^2 + (Zmesh-centerDepth).^2 < (radiusDisk-3e-3).^2;
back = Xmesh.^2 + (Zmesh-centerDepth).^2 > (radiusDisk+3e-3).^2;
fprintf('B/A inc: %.2f +/- %.2f\n', mean(estBAinst(inc),'omitnan'), ...
std(estBAinst(inc), [] ,'omitnan'));
fprintf('B/A back: %.2f +/- %.2f\n', mean(estBAinst(back),'omitnan'), ...
std(estBAinst(back), [] ,'omitnan'));

%% ADMM
% Hyperparameters BEST
tol = 1e-3;
muAlpha = 1e0; muBeta = 1e-1;
rho = 10;
maxIte = 200;

% Hyperparameters TEST
% [m,n,p] = size(bzf);
% tol = 1e-3;
% muAlpha = 1e-1; muBeta = 1e-1;
% rho = 10;
% maxIte = 200;

% Initialization
% u = [alpha0*ones(n*m,1);beta0*ones(n*m,1)];
% v = [alpha0*ones(n*m,1);beta0*ones(n*m,1)];
u = u0;
v = u;
w = zeros(size(u));

% Auxiliar
dzP = zP(2)-zP(1);
izP = round(zP./dzP);
DP = [-izP izP];
DP = spdiags(DP, [-1 0],m,m);
DP(1,:) = DP(2,:);
DP = kron(speye(n),DP);
mask = ones(n*m,1);
I = speye(2*m*n);
Ii = I(:,1:m*n);
Id = I(:,1+m*n:end);

% 
% Optimizes F(u) + R(v)

% Objective functions
Fid = []; Reg = []; Dual = [];
Fid(1) = 1/2*norm( modelFreq(u,zP,freqVec) - bzf(:) )^2;
Reg(1) = muAlpha*TVcalc_isotropic(DP*u(1:m*n),m,n,mask) + ...
    muBeta*TVcalc_isotropic(DP*u(m*n+1:end),m,n,mask);
Dual(1) = 0;
ite  = 0;
error = 1;
tic
while abs(error) > tol && ite < maxIte
    ite = ite + 1;

    % Fidelity step
    u = optimNonLinearGNFreq(zP,freqVec,bzf(:), speye(2*m*n),v-w, rho, tol,u);

    % Regularization step
    v = optimAdmmTvDy(Ii,Id,u+w, muAlpha/rho,muBeta/rho ,m,n,tol,mask,DP);
    
    % Dual update
    w = w + u - v;

    % Loss
    Fid(ite+1) = 1/2*norm( modelFreq(u,zP,freqVec) - bzf(:) )^2;
    Reg(ite+1) = muAlpha*TVcalc_isotropic(DP*v(1:m*n),m,n,mask) + ...
        muBeta*TVcalc_isotropic(DP*v(m*n+1:end),m,n,mask);
    Dual(ite+1) = norm(u-v);
    fprintf("Ite: %01i, Fid: %.3f, Reg: %.3f, Dual: %.3f\n",...
        ite,Fid(ite+1),Reg(ite+1),Dual(ite+1))
    error = Fid(ite+1) + Reg(ite+1) - Fid(ite) - Reg(ite);

    if mod(ite,10) ==1 && ite<50
        alphaArr = reshape(u(1:m*n),[m,n]);
        betaArr = reshape(u(1+m*n:end),[m,n]);
        estACtv = alphaArr*NptodB/100;
        estBAtv = 2*(betaArr-1);

        figure; tiledlayout(1,2)
        t1 = nexttile;
        imagesc(xP*1e3,zP*1e3,estACtv); colorbar;
        clim(attRange);
        title('\alpha_0 in \alpha(f) = \alpha_0 \times f^2 dB/cm')
        axis image
        colorbar;
        xlabel('Lateral distance (mm)');
        ylabel('Depth (mm)');
        hold on
        rectangle('Position',[0-radiusDisk,centerDepth-radiusDisk,...
            2*radiusDisk,2*radiusDisk]*1000, 'Curvature',1,...
            'EdgeColor','b', 'LineStyle','--', 'LineWidth',2)
        hold off

        nexttile; imagesc(xP*1e3,zP*1e3,estBAtv); colorbar;
        clim(baRange);
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

        colormap(t1,turbo); 
        pause(0.1)

    end
end
toc
%%
figure('Units','centimeters', 'Position',[5 5 10 15]),
tiledlayout(2,1)
nexttile,
plot(Fid+Reg, 'LineWidth',2)
ylabel('Objective')
xlim([3 length(Fid)])
nexttile,
plot(Dual, 'LineWidth',2)
xlabel('Number of iterations')
ylabel('Dual step norm')
xlim([3 length(Fid)])


alphaArr = reshape(DP*u(1:m*n),[m,n]);
betaArr = reshape(DP*u(1+m*n:end),[m,n]);

estACtv = alphaArr*NptodB/100;
estBAtv = 2*(betaArr-1);


[Xmesh,Zmesh] = meshgrid(xP,zP);
inc = Xmesh.^2 + (Zmesh-centerDepth).^2 < (radiusDisk-1e-3).^2;
back = Xmesh.^2 + (Zmesh-centerDepth).^2 > (radiusDisk+1e-3).^2;

fprintf('AC: %.3f +/- %.3f\n', mean(estACtv(:), 'omitnan'), ...
std(estACtv(:), [] ,'omitnan'));
fprintf('B/A: %.2f +/- %.2f\n', mean(estBAtv(:),'omitnan'), ...
std(estBAtv(:), [] ,'omitnan'));

fprintf('B/A inc: %.2f +/- %.2f\n', mean(estBAtv(inc),'omitnan'), ...
std(estBAtv(:), [] ,'omitnan'));
fprintf('B/A back: %.2f +/- %.2f\n', mean(estBAtv(back),'omitnan'), ...
std(estBAtv(:), [] ,'omitnan'));


figure; imagesc(xP*1e3,zP*1e3,estACtv); colorbar;
clim(attRange);
title('\alpha_0 in \alpha(f) = \alpha_0 \times f^2 dB/cm')
axis image
colormap turbo; colorbar;
xlabel('Lateral distance (mm)');
ylabel('Depth (mm)');
hold on
rectangle('Position',[0-radiusDisk,centerDepth-radiusDisk,...
2*radiusDisk,2*radiusDisk]*1000, 'Curvature',1,...
'EdgeColor','b', 'LineStyle','--', 'LineWidth',2)
hold off


figure; imagesc(xP*1e3,zP*1e3,estBAtv); colorbar;
clim(baRange);
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


save_all_figures_to_directory(fullfile(baseDir,'24-09-19'),'fig')
close all
%% Utility functions
function [bz,xP,zP] = getMeasurements(medium,filterParams,blockParams)
%   Gets a 2D map of measurements (mz in the IUS 2024 paper)
%   given the rf data and some hyperparameters. Downsamples in the axial
%   and lateral directions

z = medium.z;
x = medium.x;
fs = medium.fs;
rfL = medium.rfL;
rfH = medium.rfH;
rfLR = medium.rfLR;
rfHR = medium.rfHR;
v = medium.v; 

% Filtering
dz = z(2) - z(1);
freqC = filterParams.freqC;
freqTol = filterParams.freqTol;
wl = 1540/freqC/1e6;
order = round(wl*filterParams.nCycles/dz);
PLfull = getFilteredPressure(rfL,fs,freqC*1e6,freqTol*1e6,order);
PHfull = getFilteredPressure(rfH,fs,freqC*1e6,freqTol*1e6,order);
PLRfull = getFilteredPressure(rfLR,fs,freqC*1e6,freqTol*1e6,order);
PHRfull = getFilteredPressure(rfHR,fs,freqC*1e6,freqTol*1e6,order);

% Block averaging
[meanP,xP,zP] = getMeanBlock(cat(3,PLfull,PHfull,PLRfull,PHRfull),x,z,...
    blockParams);
PL = meanP(:,:,1);
PH = meanP(:,:,2);
PLR = meanP(:,:,3);
PHR = meanP(:,:,4);

% Getting measurements
attRef = medium.alphaR*freqC^2;
bz = medium.betaR*sqrt( abs(v*PL-PH)./abs(v*PLR-PHR) .*PLR./PL ).*...
    ( 1-exp(-2*attRef*zP) )/attRef./zP;

end
% ---------------------------------------------------------------------- %
% ---------------------------------------------------------------------- %
function  P = getFilteredPressure(rf,fs,freqC,freqTol,order)
%   Filters the rf data around a center frequency and returns the envelope
%   All frequencies should be in MHz

t = (0:size(rf,1)-1)'./fs;
rfMod = rf.*exp(1j*2*pi*freqC*t);

freqNyq = fs/2;
d = floor(order/2);
b = fir1(order,freqTol/freqNyq);
nCol = size(rf,2);
rfFilt = filter(b,1,[rfMod;zeros(d,nCol)]);
rfFilt = rfFilt(d+1:end,:);
P = abs(rfFilt);
end