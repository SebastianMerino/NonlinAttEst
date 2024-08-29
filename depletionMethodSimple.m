clear; close all; clc;
addpath(genpath(pwd))
baseDir = ['C:\Users\sebas\Documents\MATLAB\totalvarsimul_AC_BA\' ...
    'BA_AC_joint\rfdata'];
fileSam = 'rf_fnum3_PWNE_samBA9_att0p10f2_nc10_400kPa';
% fileSam = 'rf_fnum3_PWNE_samBA9_att0p18f2_nc10_400kPa';
fileRef = 'rf_fnum3_PWNE_refBA6_att10f2_nc10_400kPa';

% baseDir = 'C:\Users\smerino.C084288\Documents\MATLAB\BA_AC_joint\rfdata';
% fileSam = 'rf_baBack6_baInc9_att0p1.mat';
% fileSam = 'rf_ba9_attBack0p1_attInc0p18.mat';
% fileSam = 'rf_baBack6_baInc9_attBack0p1_attInc0p18.mat';
% fileRef = 'rf_ba8_att0p12_ref.mat';
NptodB = 20*log10(exp(1));

% Hyperparameters
zIni = 0.3; zFin = 5.5;
freqC = 5; freqTol = 2;
blockSize = 20; overlap = 0.8;
radiusDisk = (9)*1e-3;
centerDepth = 22.5e-3;

% Known variables
betaR = 1 + 6/2;
alphaR = 0.10*freqC.^2/NptodB*100;
alphaS = 0.1*freqC.^2/NptodB*100;
% betaR = 1 + 8/2;
% alphaR = 0.12*freqC.^2/NptodB*100;
% alphaS = 0.1*freqC.^2/NptodB*100;
v = 5; % scaling factor

%% Sample
% Loading and cropping
load(fullfile(baseDir,fileSam))
idz = z>zIni/100 & z<zFin/100;
z = z(idz);
rf1 = rf1(idz,:,1);
rf2 = rf2(idz,:,1);

% Filtering
rf1Filt = filterRfSignal(rf1,fs,freqC,freqTol);
rf2Filt = filterRfSignal(rf2,fs,freqC,freqTol);

% Getting pressures
wl = 1540/freqC/1e6;
[PL,PH,xP,zP] = getPressures(rf1Filt,rf2Filt,x,z, ...
    [blockSize*wl,blockSize*wl],overlap);

% Plotting Bmode
Bmode = db(hilbert(rf1Filt));
Bmode = Bmode - max(Bmode(:));
figure, imagesc(x*100,z*100,Bmode, [-50 0])
hold on
rectangle('Position',[0-radiusDisk,centerDepth-radiusDisk,...
    2*radiusDisk,2*radiusDisk]*100, 'Curvature',1,...
    'EdgeColor','b', 'LineStyle','--', 'LineWidth',2)
hold off
axis image
colormap gray
xlabel('Lateral [cm]')
ylabel('Depth [cm]')

%% Reference
% Loading and cropping
load(fullfile(baseDir,fileRef))
z = z(idz);
rf1 = rf1(idz,:,1);
rf2 = rf2(idz,:,1);

% Filtering
rf1Filt = filterRfSignal(rf1,fs,freqC,freqTol);
rf2Filt = filterRfSignal(rf2,fs,freqC,freqTol);

% Getting pressures
[PLR,PHR,~,~] = getPressures(rf1Filt,rf2Filt,x,z, ...
    [blockSize*wl,blockSize*wl],overlap);

PLR = mean(PLR,2);
PHR = mean(PHR,2);

%% Getting B/A
betaS = betaR*sqrt( (v*PL-PH)./(v*PLR-PHR) .*PLR./PL ).*...
    ( 1-exp(-2*alphaR*zP) )./( 1-exp(-2*alphaS*zP) ) *alphaS/alphaR;

BA = 2*(betaS-1);
cm = 100;
figure,
imagesc(xP*cm,zP*cm,BA, [5 14])
title(sprintf("B/A = %.2f \\pm %.2f",mean(BA(:)),std(BA(:))))
xlabel('Lateral [cm]')
ylabel('Depth [cm]')
axis image
colorbar
colormap pink
hold on
rectangle('Position',[0-radiusDisk,centerDepth-radiusDisk,...
    2*radiusDisk,2*radiusDisk]*100, 'Curvature',1,...
    'EdgeColor','b', 'LineStyle','--', 'LineWidth',2)
hold off

%%
function [P1,P2,xP,zP] = getPressures(rf1,rf2,x,z,blockSize,overlap)
dx = x(2) - x(1);
dz = z(2) - z(1);

% Lateral samples
wx = round(blockSize(1)*(1-overlap)/dx);  % Between windows
nx = round(blockSize(1)/dx);                 % Window size
x0 = 1:wx:length(x)-nx;
xP = x(1,x0+round(nx/2));
n  = length(x0);

% Axial samples
wz = round(blockSize(2)*(1-overlap)/dz); % Between windows
nz = round(blockSize(2)/dz); % Window size
z0 = 1:wz:length(z)-nz;
zP = z(z0 + round(nz/2))';
m  = length(z0);

P1raw = abs(hilbert(rf1));
P2raw = abs(hilbert(rf2));
P1 = zeros(m,n);
P2 = zeros(m,n);
for ii = 1:m
    for jj = 1:n
        blockP = P1raw(z0(ii):z0(ii)+nz,x0(jj):x0(jj)+nx);
        P1(ii,jj) = mean(blockP(:));

        blockP = P2raw(z0(ii):z0(ii)+nz,x0(jj):x0(jj)+nx);
        P2(ii,jj) = mean(blockP(:));
    end
end

end

function rfFilt = filterRfSignal(rf,fs,freqC,freqTol)
order = 200;
freqNyq = fs/2/1e6;
d = floor(order/2);
b = fir1(order,[freqC-freqTol freqC+freqTol]/freqNyq,'bandpass');
nCol = size(rf,2);
rfFilt = filter(b,1,[rf;zeros(d,nCol)]);
rfFilt = rfFilt(d+1:end,:);
end