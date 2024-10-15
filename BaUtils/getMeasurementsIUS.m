function [bz,zbz,xB,zB] = getMeasurementsIUS(medium,filterParams,blockParams)
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

% In case there is an additional dimension
PLfull = mean(PLfull,3);
PHfull = mean(PHfull,3);
PLRfull = mean(PLRfull,3);
PHRfull = mean(PHRfull,3);

% Downsampling lateral (moving average)
[meanP,xB] = averageLines(cat(3,PLfull,PHfull,PLRfull,PHRfull),x,blockParams);
PL = meanP(:,:,1);
PH = meanP(:,:,2);
PLR = meanP(:,:,3);
PHR = meanP(:,:,4);

% Getting measurements
if isfield(medium,'alphaR')
    attRef = medium.alphaR*(freqC/1e6)^2;
else
    attRef = medium.alphaRcoeff*(freqC/1e6)^medium.alphaRpower;
end
bz = medium.betaR*sqrt( abs(v*PL-PH)./abs(v*PLR-PHR) .*PLR./PL ).*...
    ( 1-exp(-2*attRef*z) )/attRef./z;

% Ordering in blocks
[bz, zbz] = getZBlock(bz,z,blockParams);
zB = mean(zbz,[2,3]);
end

function [bzOrd, zbz] = getZBlock(bz,z,blockParams)

% Downsampling
D = blockParams.downFactor;
z = z(1:D:end);
bz = bz(1:D:end,:);

% Parameters
dz = z(2) - z(1);
blockSize = blockParams.blockSize;
overlap = blockParams.overlap;

% Cropping
idz = z>blockParams.zlim(1) & z<blockParams.zlim(end);
z = z(idz);
bz = bz(idz,:);

% Axial samples
wz = round(blockSize(2)*(1-overlap)/dz); % Between windows
nz = round(blockSize(2)/dz); % Window size
z0 = 1:wz:length(z)-nz;
m  = length(z0);

n = size(bz,2);
bzOrd = zeros(m,n,nz);
zbz = zeros(m,n,nz);
for ii = 1:m
    for jj = 1:n
        bzBlock = bz(z0(ii):z0(ii)+nz-1,jj);
        bzOrd(ii,jj,:) = bzBlock;
        zbz(ii,jj,:) = z(z0(ii):z0(ii)+nz-1);
    end
end


end

function [Pblock,xBlock] = averageLines(P,x,blockParams)
dx = x(2) - x(1);
blockSize = blockParams.blockSize;
overlap = blockParams.overlap;

% Cropping
idx = x>blockParams.xlim(1) & x<blockParams.xlim(end);
x = x(idx); 
P = P(:,idx,:);

% Lateral samples
wx = round(blockSize(1)*(1-overlap)/dx);  % Between windows
nx = round(blockSize(1)/dx);                 % Window size
x0 = 1:wx:length(x)-nx;
xBlock = x(1,x0+round(nx/2));
n  = length(x0);

m = size(P,1);
Pblock = zeros(m,n,size(P,3));
for jj = 1:n
    blockP = P(:,x0(jj):x0(jj)+nx-1,:);
    Pblock(:,jj,:) = mean(blockP,2);
    xBlock(jj) = mean(x(x0(jj):x0(jj)+nx-1));
end
end
