%% B/A estimator with time domain approach
% by: Andres Coila

clear; close all; clc;
addpath(genpath(pwd))


%% Prepare files for sanmple and reference
baseDir = ['C:\Users\sebas\Documents\MATLAB\totalvarsimul_AC_BA\' ...
    'BA_AC_joint\rfdata'];
% fileSam = 'rf_fnum3_PWNE_samBA6_att0p10f2_nc10_400kPa';
fileSam = 'rf_fnum3_PWNE_samBA9_att0p18f2_nc10_400kPa';
fileRef = 'rf_fnum3_PWNE_refBA6_att10f2_nc10_400kPa';

% fileSam = 'rf_fnum3_SCOMP5MHz_nc10_0p10f2_saminc400_doubleangle720';
% fileSam = 'rf_fnum3_SCOMP5MHz_nc10_0p10f2_sam400_doubleangle720';
% fileRef = 'rf_fnum3_SCOMP5MHz_nc10_0p10f2_ref400_doubleangle720';

% fileSam = 'rf_baBack6_baInc9_att0p1.mat';
% fileSam = 'rf_ba9_attBack0p1_attInc0p18.mat';
% fileSam = 'rf_baBack6_baInc9_attBack0p1_attInc0p18.mat';
% fileRef = 'rf_ba8_att0p12_ref.mat';
c0 = 1540;
load(fullfile(baseDir,fileSam))
media.samLP = rf1(:,:,1);
media.samHP = rf2(:,:,1);
media.c0 = c0;
media.x = x;
media.z = z';

load(fullfile(baseDir,fileRef))
media.refLP = rf1(:,:,1);
media.refHP = rf2(:,:,1);

%% Assign hyperparameters
NptodB = 20*log10(exp(1));

param.ACsam = 0.18/NptodB;
param.APsam = 2;
param.ACref = 0.1/NptodB;
param.APref = 2;
param.BoAref = 6;

param.factor = 5;   % when P goes from 100 kPa to 1000 kPa
param.order_filter = 200; % Filter for time domain rf data
param.f0 = 5;   % MHz
param.size_wl = 10; % Moving average window for envelopes

param.width = 4; % Width to smooth envelopes laterally per B/A line
param.overlap_pct = 0.5;    % 80%
param.overlap = round((1-param.overlap_pct)*param.width);

param.col_last = param.width;
param.n = 0;

param.fs = fs;
param.fnumber = 3; % if plane wave Beamformingstill required

%% Estimation of B/A image

[BAimage, BAparam] = BAestimator(media, param);


figure; imagesc(BAparam.x, BAparam.z, BAimage);
font = 20;
axis image
colormap pink; colorbar;
clim([5 10])
%font = 20;
set(gca,'fontsize',font)
xlabel('\bfLateral distance (mm)');
ylabel('\bfDepth (mm)');
%title("\bf BoA map");
set(gca,'FontSize',font);
ylim([5 55])

figure; plot(BAparam.z, BAimage(:,ceil(size(BAimage,2)/2)),'LineWidth',3);
grid on;
set(gca,'fontsize',font)
xlabel('\bfDepth (mm)');
ylabel('\bf B/A');
%  title("\bf B/A");
set(gca,'FontSize',font);
ylim([0 15])

BAeff.image = BAimage; %QUS
BAeff.lateral = BAparam.x*1e-3;
BAeff.axial = BAparam.z*1e-3;
BAeff.mz = BAparam.mz;

% filename = ['FULLMAPv60_',fileSam,'_',fileRef,'.mat'];
% save(fullfile(['C:\Users\sebas\Documents\MATLAB\totalvarsimul_AC_BA\' ...
%     'BA_AC_joint\maps'],filename),'BAeff');

%% Getting system
freq = 5;
muAlpha = 1; muBeta = 1;
zlim = [5 55]*1e-3;
nwl = 40;
beta0 = 1+(10.5)/2;
alpha0 = 0.12*5^2/8.686*100; % dB/cm -> Np/m

mzaxis = BAeff.mz;
x = BAeff.lateral';
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
maxIte = 400;
muAlpha = 10;
muBeta = 10;
% muAlpha = 0; muBeta = 0;
beta0 = 1+(8.5)/2;
alpha0 = 0.1*freq^2; % dB/cm -> Np/m

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
xlim([5 length(loss)])
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
pause(0.5)




%% Getting mz
function [BAimage, BAparam] = BAestimator(media, param)
scale = param.factor;

ACref = param.ACref;
ACsam = param.ACsam;
alpha_power_sam = param.APsam;
alpha_power_ref = param.APref;

c0 = media.c0;
x = media.x;
z = media.z;

%% Filtering
samLPfull = media.samLP;
samHPfull = media.samHP;
refLPfull = media.refLP;
refHPfull = media.refHP;

filterParams.order = param.order_filter;
filterParams.fs = param.fs;
filterParams.freqC = param.f0;
filterParams.freqTol = 1.5;
samLPfull = filterRfSignal(samLPfull,filterParams);
samHPfull = filterRfSignal(samHPfull,filterParams);
refLPfull = filterRfSignal(refLPfull,filterParams);
refHPfull = filterRfSignal(refHPfull,filterParams);

%%

f0 = param.f0;
dz = z(2)-z(1);
L = param.size_wl * (c0/(f0*1e6)) / dz;
f0_fund_sam = param.f0;
f0_fund_ref = param.f0;

%% BoA estimation

envLPsam = abs(hilbert(samLPfull));
envHPsam = abs(hilbert(samHPfull));
envLPref = abs(hilbert(refLPfull));
envHPref = abs(hilbert(refHPfull));

envLPsam = movmean(envLPsam,[0 param.width-1],2);
envHPsam = movmean(envHPsam,[0 param.width-1],2);

envLPsam = envLPsam(:,1:param.overlap:end);
envHPsam = envHPsam(:,1:param.overlap:end);
xBA = x(param.overlap:param.overlap:end);

envLPref = mean(envLPref,2);
envHPref = mean(envHPref,2);

envLPsam = movmean(envLPsam,L);
envHPsam = movmean(envHPsam,L);
envLPref = movmean(envLPref,L);
envHPref = movmean(envHPref,L);

GAPsam = scale*envLPsam - envHPsam;
GAPref = scale*envLPref - envHPref;
% GAPsam = movmean(GAPsam,L);
% GAPref = movmean(GAPref,L);

alpha1 = ACsam;

ATT1MAPsam = ...
    1./exp( mean( alpha1.* f0_fund_sam.^alpha_power_sam .*(z*1e2) , 2 ) );
ATT1ref = 1./exp( ACref.*(f0_fund_ref.^alpha_power_ref) .*(z*1e2));

P22sam_smooth = sqrt(abs(GAPsam).*envLPref);
P22ref_smooth = sqrt(abs(GAPref).*envLPsam);

alpha1_sam =  mean( ACsam .* f0_fund_sam.^alpha_power_sam , 2) *100;
alpha1_ref =   ( ACref.* ( (f0_fund_ref).^alpha_power_ref)  )*100; % Np/m

mz = ((1+0.5*param.BoAref))* (P22sam_smooth./P22ref_smooth) ...
    .* (1- ATT1ref.*ATT1ref)./(alpha1_ref.*(z));

ratio_beta = (P22sam_smooth./P22ref_smooth) ...
    .* (alpha1_sam.*z./(alpha1_ref.*z)) ...
    .* (1- ATT1ref.*ATT1ref)./(1- ATT1MAPsam.*ATT1MAPsam);

BoAsam = ( ratio_beta*(1+0.5*param.BoAref) - 1 )*2;


BoAsam_matrix = BoAsam;
mzmatrix = mz;

BAimage = BoAsam_matrix;
BAparam.mz =  mzmatrix;
BAparam.z = z(:)*1e3;
BAparam.x = xBA(:)*1e3;

end

%% Utility functions
function rfFilt = filterRfSignal(rf,filterParams)
% fs,freqC,freqTol
order = filterParams.order;
freqNyq = filterParams.fs/2/1e6;
fL = filterParams.freqC - filterParams.freqTol;
fH = filterParams.freqC + filterParams.freqTol;

d = floor(order/2);
b = fir1(order,[fL fH]/freqNyq,'bandpass');
nCol = size(rf,2);
rfFilt = filter(b,1,[rf;zeros(d,nCol)]);
rfFilt = rfFilt(d+1:end,:);
end
