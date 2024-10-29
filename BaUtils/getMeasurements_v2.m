function [bz,xP,zP] = getMeasurements_v2(medium,filterParams,blockParams,varargin)
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
DEPfull = getFilteredPressure(v*rfL-rfH,fs,freqC,freqTol,order);
PLRfull = getFilteredPressure(rfLR,fs,freqC,freqTol,order);
DEPRfull = getFilteredPressure(v*rfLR-rfHR,fs,freqC,freqTol,order);

% In case there is an additional dimension
PLfull = mean(PLfull,3);
DEPfull = mean(DEPfull,3);
PLRfull = mean(PLRfull,3);
DEPRfull = mean(DEPRfull,3);

% Block averaging
[meanP,xP,zP] = getMeanBlock(cat(3,PLfull,DEPfull,PLRfull,DEPRfull),x,z,...
    blockParams);
PL = meanP(:,:,1);
DEP = meanP(:,:,2);
PLR = meanP(:,:,3);
DEPR = meanP(:,:,4);

% Getting measurements
if nargin==3
    attRef = medium.alphaR*(freqC/1e6)^2;
else
    attRef = medium.alphaR*(freqC/1e6).^varargin{1};
end
bz = medium.betaR*sqrt( abs(DEP)./abs(DEPR) .*PLR./PL ).*...
    ( 1-exp(-2*attRef*zP) )/attRef./zP;
end
