% My implementation of the new method with frequency compounding
% Used for testing several frequencies

clear; close all; clc;
addpath(genpath(pwd))
baseDir = ['C:\Users\sebas\Documents\MATLAB\totalvarsimul_AC_BA\' ...
    'BA_AC_joint\rfdata'];
% baseDir = 'C:\Users\smerino.C084288\Documents\MATLAB\BA_AC_joint\rfdata';

% fileSam = 'rf_fnum3_PWNE_samBA9_att0p10f2_nc10_400kPa';
% fileSam = 'rf_fnum3_PWNE_samBA9_att0p18f2_nc10_400kPa';
% fileRef = 'rf_fnum3_PWNE_refBA6_att10f2_nc10_400kPa';

% fileSam = 'rf_fnum3_SCOMP5MHz_nc10_0p10f2_saminc400_doubleangle720';
% fileSam = 'rf_fnum3_SCOMP5MHz_nc10_0p10f2_sam400_doubleangle720';
% fileRef = 'rf_fnum3_SCOMP5MHz_nc10_0p10f2_ref400_doubleangle720';

% fileSam = 'rf_baBack6_baInc9_att0p1.mat';
% fileSam = 'rf_ba9_attBack0p1_attInc0p18.mat';
% fileSam = 'rf_baBack6_baInc9_attBack0p1_attInc0p18.mat';
% fileRef = 'rf_ba8_att0p12_ref.mat';

fileSam = 'RFfn3_PWNE_samBA9_att0p10f2_nc10_400kPa';
fileRef = 'RFfn3_PWNE_refBA6_att12f2_nc10_400kPa';

% Auxiliar variables
NptodB = 20*log10(exp(1));
radiusDisk = (9)*1e-3;
centerDepth = 22.5e-3;

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
blockParams.blockSize = [10 20]*wl; 
blockParams.overlap = 0.8;
blockParams.zlim = [0.3; 5.5]/100;
blockParams.xlim = [-2; 2]/100;

%% Getting B/A
bzf = [];
freqVec = [4,4.5,5,5.5,6];
for iFreq = 1:length(freqVec)
    filterParams.freqC = freqVec(iFreq);
    [bz,xP,zP] = getMeasurements(medium,filterParams,blockParams);
    bzf(:,:,iFreq) = bz;
end
%% Measurements
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
betaIni = 1+(7.5)/2;
alphaIni = 0.18/NptodB*100;

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
%%
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

figure; imagesc(xP*1e3,zP*1e3,estBAinst); colorbar;
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

%%
function [bz,xP,zP] = getMeasurements(medium,filterParams,blockParams)

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
% disp(wl*filterParams.nCycles)
order = round(wl*filterParams.nCycles/dz);
PLfull = getFilteredPressure(rfL,fs,freqC*1e6,freqTol*1e6,order);
PHfull = getFilteredPressure(rfH,fs,freqC*1e6,freqTol*1e6,order);
PLRfull = getFilteredPressure(rfLR,fs,freqC*1e6,freqTol*1e6,order);
PHRfull = getFilteredPressure(rfHR,fs,freqC*1e6,freqTol*1e6,order);

% figure,
% tiledlayout(1,4)
% nexttile,
% imagesc(x*100,z*100,PLfull)
% title("Measurements b(z,f)")
% xlabel('Lateral [cm]')
% ylabel('Depth [cm]')
% axis image
% colorbar
% colormap parula
% 
% nexttile,
% imagesc(x*100,z*100,PHfull)
% title("Measurements b(z,f)")
% xlabel('Lateral [cm]')
% ylabel('Depth [cm]')
% axis image
% colorbar
% colormap parula
% 
% nexttile,
% imagesc(x*100,z*100,PLRfull)
% title("Measurements b(z,f)")
% xlabel('Lateral [cm]')
% ylabel('Depth [cm]')
% axis image
% colorbar
% colormap parula
% 
% nexttile,
% imagesc(x*100,z*100,PHRfull)
% title("Measurements b(z,f)")
% xlabel('Lateral [cm]')
% ylabel('Depth [cm]')
% axis image
% colorbar
% colormap parula


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