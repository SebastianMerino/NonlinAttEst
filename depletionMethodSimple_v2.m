% My implementation of the depletion method with my own functions
% Used for testing several frequencies

clear; close all; clc;
addpath(genpath(pwd))
baseDir = ['C:\Users\sebas\Documents\MATLAB\totalvarsimul_AC_BA\' ...
    'BA_AC_joint\rfdata'];
% baseDir = 'C:\Users\smerino.C084288\Documents\MATLAB\BA_AC_joint\rfdata';

fileSam = 'rf_fnum3_PWNE_samBA9_att0p10f2_nc10_400kPa';
% fileSam = 'rf_fnum3_PWNE_samBA9_att0p18f2_nc10_400kPa';
fileRef = 'rf_fnum3_PWNE_refBA6_att10f2_nc10_400kPa';

% fileSam = 'rf_fnum3_SCOMP5MHz_nc10_0p10f2_saminc400_doubleangle720';
% fileSam = 'rf_fnum3_SCOMP5MHz_nc10_0p10f2_sam400_doubleangle720';
% fileRef = 'rf_fnum3_SCOMP5MHz_nc10_0p10f2_ref400_doubleangle720';

% fileSam = 'rf_baBack6_baInc9_att0p1.mat';
% fileSam = 'rf_ba9_attBack0p1_attInc0p18.mat';
% fileSam = 'rf_baBack6_baInc9_attBack0p1_attInc0p18.mat';
% fileRef = 'rf_ba8_att0p12_ref.mat';

% Auxiliar variables
NptodB = 20*log10(exp(1));
radiusDisk = (9)*1e-3;
centerDepth = 22.5e-3;

% Hyperparameters
freqC = 5; 
zIni = 0.3; zFin = 5.5;
wl = 1540/freqC/1e6;
v = 5; % scaling factor

% Known variables
betaR = 1 + 6/2;
alphaR = 0.10*freqC.^2/NptodB*100;
alphaS = 0.10*freqC.^2/NptodB*100;
% betaR = 1 + 8/2;
% alphaR = 0.12*freqC.^2/NptodB*100;
% alphaS = 0.1*freqC.^2/NptodB*100;

%% Loading and cropping rf data
% Sample
sample = load(fullfile(baseDir,fileSam));
idz = sample.z>zIni/100 & sample.z<zFin/100;
z = sample.z(idz)';
x = sample.x;
fs = sample.fs;
rfL = sample.rf1(idz,:,1);
rfH = sample.rf2(idz,:,1);
% Reference
ref = load(fullfile(baseDir,fileRef));
rfLR = ref.rf1(idz,:,1);
rfHR = ref.rf2(idz,:,1);


%% Filtering
freqTol = 1;
nCycles = 10; % Number of cycles of the initial filter

dz = z(2) - z(1);
order = round(wl*nCycles/dz);
PLfull = getFilteredPressure(rfL,fs,freqC*1e6,freqTol*1e6,order);
PHfull = getFilteredPressure(rfH,fs,freqC*1e6,freqTol*1e6,order);
PLRfull = getFilteredPressure(rfLR,fs,freqC*1e6,freqTol*1e6,order);
PHRfull = getFilteredPressure(rfHR,fs,freqC*1e6,freqTol*1e6,order);


% Plotting Bmode
Bmode = db(PLfull);
Bmode = Bmode - max(Bmode(:));
figure, imagesc(x*100,z*100,Bmode, [-50 0])
% hold on
% rectangle('Position',[0-radiusDisk,centerDepth-radiusDisk,...
%     2*radiusDisk,2*radiusDisk]*100, 'Curvature',1,...
%     'EdgeColor','b', 'LineStyle','--', 'LineWidth',2)
% hold off
axis image
colormap gray
xlabel('Lateral [cm]')
ylabel('Depth [cm]')

%% 
blockSize = 20; overlap = 0.8;
[meanP,xP,zP] = getMeanBlock(cat(3,PLfull,PHfull,PLRfull,PHRfull),x,z,...
    [blockSize*wl,blockSize*wl],overlap);
PL = meanP(:,:,1);
PH = meanP(:,:,2);
PLR = meanP(:,:,3);
PHR = meanP(:,:,4);

figure,imagesc(x*100,z*100,PL)
axis image
colorbar
%% Getting B/A
betaS = betaR*sqrt( abs(v*PL-PH)./abs(v*PLR-PHR) .*PLR./PL ).*...
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
% hold on
% rectangle('Position',[0-radiusDisk,centerDepth-radiusDisk,...
%     2*radiusDisk,2*radiusDisk]*100, 'Curvature',1,...
%     'EdgeColor','b', 'LineStyle','--', 'LineWidth',2)
% hold off

%% 
mz = betaR*sqrt( abs(v*PL-PH)./abs(v*PLR-PHR) .*PLR./PL ).*...
    ( 1-exp(-2*alphaR*zP) )./alphaR./zP;
figure,
imagesc(xP*cm,zP*cm,mz)
xlabel('Lateral [cm]')
ylabel('Depth [cm]')
axis image
colorbar
colormap parula
%%
function  P = getFilteredPressure(rf,fs,freqC,freqTol,order)
t = (0:size(rf,1)-1)'./fs;
rfMod = rf.*exp(1j*2*pi*freqC*t);

% figure, plot(mean(abs(fft(rf)),2))
% figure, plot(mean(abs(fft(rfMod)),2))

freqNyq = fs/2;
d = floor(order/2);
b = fir1(order,freqTol/freqNyq);
nCol = size(rf,2);
rfFilt = filter(b,1,[rfMod;zeros(d,nCol)]);
rfFilt = rfFilt(d+1:end,:);
P = abs(rfFilt);
end


function [Pblock,xBlock,zBlock] = getMeanBlock(P,x,z,blockSize,overlap)
dx = x(2) - x(1);
dz = z(2) - z(1);

% Lateral samples
wx = round(blockSize(1)*(1-overlap)/dx);  % Between windows
nx = round(blockSize(1)/dx);                 % Window size
x0 = 1:wx:length(x)-nx;
xBlock = x(1,x0+round(nx/2));
n  = length(x0);

% Axial samples
wz = round(blockSize(2)*(1-overlap)/dz); % Between windows
nz = round(blockSize(2)/dz); % Window size
z0 = 1:wz:length(z)-nz;
zBlock = z(z0 + round(nz/2));
m  = length(z0);

Pblock = zeros(m,n,size(P,3));
for ii = 1:m
    for jj = 1:n
        blockP = P(z0(ii):z0(ii)+nz,x0(jj):x0(jj)+nx,:);
        Pblock(ii,jj,:) = mean(blockP,[1 2]);
    end
end

end