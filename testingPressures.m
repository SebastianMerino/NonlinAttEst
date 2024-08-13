clear; close all; clc;
addpath(genpath(pwd))

fileSam = 'rf_fnum3_PWNE_samBA12_att0p18f2_nc10_400kPa';
fileRef = 'rf_fnum3_PWNE_refBA6_att10f2_nc10_400kPa';

load(fileSam)
idz = z>0.16/100 & z<5.7/100;
z = z(idz);
rf1 = rf1(idz,:,1);
rf2 = rf2(idz,:,1);

cm = 100;
figure, imagesc(x*cm,z*cm,db(hilbert(rf1)))
colormap gray
axis image

%% Spectrum of signals
m = size(rf1,1);
Sp1 = mean(abs(fft(rf1)),2);
Sp2 = mean(abs(fft(rf2)),2);
f = (0:m-1)/m*fs/1e6;
figure,
plot(f,db([Sp1,Sp2]))
xlim([0 20])
legend('rf1', 'rf2')
%% Filtering
freqC = 5;
freqTol = 2;
rf1Filt = filterRfSignal(rf1,fs,freqC,freqTol);
rf2Filt = filterRfSignal(rf2,fs,freqC,freqTol);

Sp1Filt = mean(abs(fft(rf1Filt)),2);
plot(f,db([Sp1,Sp1Filt]))
xlim([0 20])
legend('Original', 'Filtered')

%% Getting pressure
blockSize = 30;
overlap = 0.8;

wl = 1540/freqC/1e6;
dx = x(2) - x(1);
dz = z(2) - z(1);

% Lateral samples
wx = round(blockSize*wl*(1-overlap)/dx);  % Between windows
nx = round(blockSize*wl/dx);                 % Window size
x0 = 1:wx:length(x)-nx;
xBlock = x(1,x0+round(nx/2));
n  = length(x0);

% Axial samples
wz = round(blockSize*wl*(1-overlap)/dz); % Between windows
nz = round(blockSize*wl/dz); % Window size
z0 = 1:wz:length(z)-nz;
zBlock = z(z0 + round(nz/2));
m  = length(z0);

P1raw = abs(hilbert(rf1Filt));
P2raw = abs(hilbert(rf2Filt));
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
%% Testing functions
[PL,PH,xP,zP] = getPressures(fileSam);
[PLR,PHR,~,~] = getPressures(fileRef);

NptodB = 20*log10(exp(1));
betaR = 1 + 6/2;
alphaR = 0.1*freqC.^2/NptodB*100;
alphaS = 0.18*freqC.^2/NptodB*100;
v = 5; % scaling factor

betaS = betaR*sqrt( (v*PL-PH)./(v*PLR-PHR) .*PLR./PL ).*...
    ( 1-exp(-2*alphaR*zP) )./( 1-exp(-2*alphaS*zP) ) *alphaS/alphaR;

BA = 2*(betaS-1);
figure,
imagesc(xP*cm,zP*cm,BA, [6 14])
axis image
colorbar
colormap pink
%%
function [P1,P2,xP,zP] = getPressures(filename)

% Hyperparameters
zIni = 0.2; zFin = 5.5;
freqC = 5;
freqTol = 2;
blockSize = 40;
overlap = 0.8;

load(filename, 'rf1','rf2','x','z','fs')
idz = z>zIni/100 & z<zFin/100;
z = z(idz);
rf1 = rf1(idz,:,1);
rf2 = rf2(idz,:,1);

% Filtering
rf1Filt = filterRfSignal(rf1,fs,freqC,freqTol);
rf2Filt = filterRfSignal(rf2,fs,freqC,freqTol);

% Getting pressure
wl = 1540/freqC/1e6;
dx = x(2) - x(1);
dz = z(2) - z(1);

% Lateral samples
wx = round(blockSize*wl*(1-overlap)/dx);  % Between windows
nx = round(blockSize*wl/dx);                 % Window size
x0 = 1:wx:length(x)-nx;
xP = x(1,x0+round(nx/2));
n  = length(x0);

% Axial samples
wz = round(blockSize*wl*(1-overlap)/dz); % Between windows
nz = round(blockSize*wl/dz); % Window size
z0 = 1:wz:length(z)-nz;
zP = z(z0 + round(nz/2))';
m  = length(z0);

P1raw = abs(hilbert(rf1Filt));
P2raw = abs(hilbert(rf2Filt));
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