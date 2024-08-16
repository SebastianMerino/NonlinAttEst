clear; close all; clc;
addpath(genpath(pwd))
% baseDir = ['C:\Users\sebas\Documents\MATLAB\totalvarsimul_AC_BA\' ...
%     'BA_AC_joint\rfdata'];
fileSam = 'rf_fnum3_PWNE_samBA12_att0p18f2_nc10_400kPa';
fileRef = 'rf_fnum3_PWNE_refBA6_att10f2_nc10_400kPa';

baseDir = 'C:\Users\smerino.C084288\Documents\MATLAB\BA_AC_joint\rfdata';
% fileSam = 'rf_baBack6_baInc9_att0p1.mat';
% fileSam = 'rf_ba9_attBack0p1_attInc0p18.mat';
% fileSam = 'rf_baBack6_baInc9_attBack0p1_attInc0p18.mat';
% fileRef = 'rf_ref_ba8_att0p12.mat';
NptodB = 20*log10(exp(1));

% Hyperparameters
zIni = 0.3; zFin = 5.5;
freqC = 5; freqTol = 2;
blockSize = 20; overlap = 0.8;
radiusDisk = (9)*1e-3;
centerDepth = 22.5e-3;

% Known variables
betaR = 1 + 6/2;
alphaR = 0.1*freqC.^2/NptodB*100;
v = 5; % scaling factor

%% Sample
% Loading and cropping
load(fullfile(baseDir,fileSam))
idz = z>zIni/100 & z<zFin/100;
z = z(idz)';
rf1 = rf1(idz,:,1);
rf2 = rf2(idz,:,1);

% Filtering
rf1Filt = filterRfSignal(rf1,fs,freqC,freqTol);
rf2Filt = filterRfSignal(rf2,fs,freqC,freqTol);

PL = abs(hilbert(rf1Filt));
PH = abs(hilbert(rf2Filt));

% Plotting Bmode
Bmode = db(hilbert(rf1Filt));
Bmode = Bmode - max(Bmode(:));
figure, imagesc(x*100,z*100,Bmode, [-50 0])
axis image
colormap gray
xlabel('Lateral [cm]')
ylabel('Depth [cm]')

%% Reference
% Loading and cropping
load(fullfile(baseDir,fileRef))
z = z(idz)';
rf1 = rf1(idz,:,1);
rf2 = rf2(idz,:,1);

% Filtering
rf1Filt = filterRfSignal(rf1,fs,freqC,freqTol);
rf2Filt = filterRfSignal(rf2,fs,freqC,freqTol);

% Getting pressures
PLR = abs(hilbert(rf1Filt));
PHR = abs(hilbert(rf2Filt));


%% Filtering pressures
pulseLength = 10;
downFactorLat = 2;

wl = 1540/freqC/1e6;
dz = z(2)-z(1);
nz = round(pulseLength*wl/dz /2)*2-1;

PL = movmean(PL, nz);
PH = movmean(PH, nz);
PLR = movmean(PLR, nz);
PHR = movmean(PHR, nz);

PL = movmean(PL, [0 1],downFactorLat);
PL = PL(:,1:downFactorLat:end);
PH = movmean(PH, [0 1],downFactorLat);
PH = PH(:,1:downFactorLat:end);
PLR = mean(PLR, 2);
PHR = mean(PHR, 2);

x = x(1:2:end);
mzaxis = betaR*( 1-exp(-2*alphaR.*z) )./alphaR.*z .*...
    sqrt( (v*PL-PH)./(v*PLR-PHR) .*PLR./PL );

%%
figure, imagesc(x*100,z*100, real(mzaxis))
axis image
%% Getting system
freq = 5;
muAlpha = 1/100; muBeta = 1;
% muAlpha = 1; muBeta = 1;
nwl = 40;
beta0 = 1+(10.5)/2;
alpha0 = 0.1*5^2/8.686*100; % dB/cm -> Np/m

dz = z(2)-z(1);
dx = x(2)-x(1);

wl = 1540/1e6/freq;
blocksize = nwl*wl; % Change to test different sizes
nz = floor(blocksize/dz);
nx = floor(blocksize/dx);

% zlim = [0 80]*1e-3;
zlim = [10 42]*1e-3;

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
nrow = length(zfin);
ncol = L2;

%% MY VERSION
% Constructing arrays
zblock = (z(zini(1)):1e-3:z(zfin(1))); % m
np = length(zblock); % Number of points per block

mzBA = zeros(np,nrow*ncol);
zBA = zeros(np,nrow*ncol); % zBA
for pixj = 1:ncol
    for pixi = 1:nrow
        zblock = (z(zini(pixi)):1e-3:z(zfin(pixi))); % m
        for zi = 1:length(zblock)
            [~, idz] = min(abs(zblock(zi)-z));
            zBA(zi,pixi+(pixj-1)*nrow) = z(idz);
            mzBA(zi,pixi+(pixj-1)*nrow) = mzaxis(idz,pixj);
            % X((pixj-1)*np+(pixi-1)*ncol*np+zi) = z(idz);
            % Y((pixj-1)*np+(pixi-1)*ncol*np+zi) = mzaxis(idz,pixj);
        end
    end
end

X = zBA(:); % zBA
Y = mzBA(:); % mzBA

% Iterating
theta = [alpha0*ones(ncol*nrow,1);beta0*ones(ncol*nrow,1)];
regMatrix = blkdiag(muAlpha*speye(ncol*nrow),muBeta*speye(ncol*nrow));
for ite =1:100
    jcb = jacobian(theta,X);
    res = Y - model(theta,X);
    [step,~] = cgs(jcb'*jcb + regMatrix,jcb'*-res);
    % step = (jcb'*jcb)\(jcb'*-res);
    theta = theta + step;
end

alphaArr = theta(1:ncol*nrow);
betaArr = theta(ncol*nrow+1:end);
alpha_dB = alphaArr/5^2/100*8.686; % Np/cm/MH^2 : 0.1

estAC2 = reshape(alpha_dB,[nrow,ncol]);
estBA2 = reshape(2*(betaArr-1),[nrow,ncol]);

disp(['AC: ',num2str(median(estAC2(:), 'omitnan'))]);
disp(['BA: ',num2str(median(estBA2(:), 'omitnan'))]);

figure; imagesc(x*1e3,z*1e3,estAC2); colorbar;
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
pause(0.5)
figure; imagesc(x*1e3,z*1e3,estBA2); colorbar;
clim([6 15]);
title('B/A');
title(['B/A = ',num2str(median(estBA2(:),'omitnan'),'%.1f')]);
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

%%
function rfFilt = filterRfSignal(rf,fs,freqC,freqTol)
order = 200;
freqNyq = fs/2/1e6;
d = floor(order/2);
b = fir1(order,[freqC-freqTol freqC+freqTol]/freqNyq,'bandpass');
nCol = size(rf,2);
rfFilt = filter(b,1,[rf;zeros(d,nCol)]);
rfFilt = rfFilt(d+1:end,:);
end

function jcb = jacobian(theta,x)
mn = size(theta,1)/2; % number of points in image
p = length(x)/mn;

alphaArr = theta(1:mn)';
betaArr = theta(mn+1:end)';
zBA = reshape(x,[p,mn]);

dMLfda = -(betaArr).*((2*zBA.*alphaArr+1).*exp(-2*alphaArr.*zBA) - 1)./(alphaArr.^2.*zBA);
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
