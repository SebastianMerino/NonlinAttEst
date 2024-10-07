% My implementation of the new method with frequency compounding
% Used for testing several frequencies

clear; close all; clc;
addpath(genpath(pwd))
baseDir = 'C:\Users\sebas\Documents\Data\Nonlinearity\homoAttMismatch\f2';
% baseDir = 'C:\Users\sebas\Documents\Data\Nonlinearity\rfdata';

% Auxiliar variables
NptodB = 20*log10(exp(1));
radiusDisk = (9)*1e-3;
centerDepth = 22.5e-3;

baRange = [5 13];
attRange = [0.05,0.2];

%% Preparing data
% Known variables
medium.v = 5;
medium.betaR = 1 + 6/2;
medium.alphaR = 0.1/NptodB*100; % alpha0 in dB/100/MHz2

% Filtering parameters
filterParams.freqC = 5;
filterParams.freqTol = 0.5;
filterParams.nCycles = 10; % Number of cycles of the initial filter

% Subsampling parameters
wl = 1540/5e6; % Mean central frequency
blockParams.blockSize = [40 40]*wl;
blockParams.overlap = 0.8;
blockParams.zlim = [0.8; 5]/100;
blockParams.xlim = [-2.5; 2.5]/100;

freqVec = [4,5,6]; % FRECUENCIES FOR FILTERING

%% Getting B/A
bzf = [];
iFreq = 1;


switch iFreq
    case 1
        fileSam = 'RFfn3_PWNE4MHz_samBA9_att0p1f2_nc10_400kPa';
        fileRef = 'RFfn3_PWNE4MHz_refBA6_att0p10f2p0_400kPa';
    case 2
        fileSam = 'RFfn3_PWNE_samBA9_att0p10f2_nc10_400kPa';
        fileRef = 'RFfn3_PWNE_refBA6_att10f2_nc10_400kPa';
    case 3
        fileSam = 'RFfn3_PWNE6MHz_samBA9_att0p1f2_nc10_400kPa';
        fileRef = 'RFfn3_PWNE6MHz_refBA6_att0p10f2p0_400kPa';
end
% Sample
sample = load(fullfile(baseDir,fileSam));
medium.z = sample.z';
medium.x = sample.x;
medium.fs = sample.fs;
medium.rfL = sample.rf1(:,:,:);
medium.rfH = sample.rf2(:,:,:);
clear sample

% Reference
ref = load(fullfile(baseDir,fileRef));
medium.rfLR = ref.rf1(:,:,:);
medium.rfHR = ref.rf2(:,:,:);
clear ref

filterParams.freqC = freqVec(iFreq);

%%
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
freqC = filterParams.freqC*1e6;
freqTol = filterParams.freqTol*1e6;
order = round(filterParams.nCycles/freqC*fs);
PLfull = getFilteredPressure(rfL,fs,freqC,freqTol,order);
PHfull = getFilteredPressure(rfH,fs,freqC,freqTol,order);
PLRfull = getFilteredPressure(rfLR,fs,freqC,freqTol,order);
PHRfull = getFilteredPressure(rfHR,fs,freqC,freqTol,order);

%% Denoising algorithms

idx = x>blockParams.xlim(1) & x<blockParams.xlim(end);
idz = z>blockParams.zlim(1) & z<blockParams.zlim(end);
xCrop = x(idx); zCrop = z(idz);
PLtest = PLfull(idz,idx,1);
PHtest = PHfull(idz,idx,1);

figure,
imagesc(xCrop*100,zCrop*100,db(PLtest))
axis image
colorbar

% PLdenoised = wdenoise2(PLfull(1:20:end,:),7);
PLdenoised = imnlmfilt(sqrt(PLtest(1:20:end,:)),'DegreeOfSmoothing',20).^2;

figure,
imagesc(xCrop*100,zCrop*100,db(PLdenoised))
axis image
colorbar

% figure,
% imagesc(x*100,z*100,PLfull-PLdenoised)
% axis image
% colorbar


%% Block averaging
% In case there is an additional dimension
PLfull = mean(PLfull,3);
PHfull = mean(PHfull,3);
PLRfull = mean(PLRfull,3);
PHRfull = mean(PHRfull,3);

% Block averaging
[meanP,xP,zP] = getMeanBlock(cat(3,PLfull,PHfull,PLRfull,PHRfull),x,z,...
    blockParams);
% meanP = meanP.^2;
PL = meanP(:,:,1);
PH = meanP(:,:,2);
PLR = meanP(:,:,3);
PHR = meanP(:,:,4);

% Getting measurements
attRef = medium.alphaR*(freqC/1e6)^2;
bz = medium.betaR*sqrt( abs(v*PL-PH)./abs(v*PLR-PHR) .*PLR./PL ).*...
    ( 1-exp(-2*attRef*zP) )/attRef;

%%
figure, surf(xP*100,zP*100,bz)
% axis image
colorbar
%%
figure, imagesc(xP*100,zP*100,diff(bz))
axis image
colorbar

bzNew = log(abs(diff(bz)));
figure, imagesc(xP*100,zP*100,bzNew)
axis image
colorbar


bzNew = medfilt2(bzNew,[3 3],'symmetric');
figure, imagesc(xP*100,zP(2:end)*100,bzNew)
axis image
colorbar


figure,plot(zP(2:end)*100,mean(bzNew,2))

% ---------------------------------------------------------------------- %
% ---------------------------------------------------------------------- %
