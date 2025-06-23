% New methods with frequency compounding. Requires three transmissions
% Used for showing results

setup;
baseDir = 'Q:\smerino\Nonlinearity\phantom';
resultsDir = 'Q:\smerino\Nonlinearity\resultsJASA\phantomDepletion';
[~,~,~] = mkdir(resultsDir);

% Auxiliar variables
NptodB = 20*log10(exp(1));

rz = 1.5e-2; rx = 1.3e-2;
cz = 4.2e-2; cx = 1.5e-2;

baRange = [4 12];
attRange = [0,1];
ylimBm = [2 5.3];
imPosition = [100 200 250 250];
alphaInit = 0.076;
baInit = 5.4;

c0 = 1470;
%% Preparing data
freqVec = [5,6,7]; % FRECUENCIES FOR FILTERING

% Known variables
medium.v = 5;

% % Reference adriana
medium.betaR = 1 + 5.4/2;
medium.alphaRcoeff = 0.3504/NptodB*100; %alpha0 in dB/m/MHz^2
medium.alphaRpower = 1.1483;

% Filtering parameters
filterParams.freqC = mean(freqVec);
filterParams.freqTol = 0.5;
filterParams.nCycles = 10; % Number of cycles of the initial filter

% Subsampling parameters
wl = c0/filterParams.freqC/1e6; % Mean central frequency
blockParams.blockSize = [30 30]*wl;
blockParams.overlap = 0.8;
blockParams.zlim = [2.4; 5]/100;
blockParams.xlim = [-1.8; 2]/100;

% Attenuation estimation parameters
blockParamsAcs.zInf = 2.4;
blockParamsAcs.zSup = 5;
blockParamsAcs.xInf = -1.8;
blockParamsAcs.xSup = 2;
blockParamsAcs.blocksize = [30 30]*wl;
blockParamsAcs.overlap = 0.8;
freqL = 4.5e6; freqH = 6.5e6;

% Regularization parameters
tol = 1e-3;
muRsld = 10^3;
muBswift = 10^3; muCswift = 10^0.5;
ratioCutOff = 10;
reject = 0.1;
extension = 3;



blockParams.downFactor = 20;

%% Measurements
medium.fs = 30e6; % new sampling frequency
medium.maxDepth = 6e-2;
dz = (1/medium.fs)*c0/2;
medium.z = (0:dz:medium.maxDepth)';

freq = 6;
fileSam = "RFfn2_L11-5v_8v40v_"+freq+"MHz_inc_cornstachx3_2";
fileRef = "RFfn2_L11-5v_8v40v_"+freq+"MHz_ref_2";

% Sample
sample = load(fullfile(baseDir,fileSam));
medium.x = sample.x;
q = 1000; p = round(q*medium.fs/sample.fs);
medium.rfL = resample(sample.rf1(:,:,:),p,q);
medium.rfH = resample(sample.rf2(:,:,:),p,q);
medium.rfL = medium.rfL(1:length(medium.z),:,:);
medium.rfH = medium.rfH(1:length(medium.z),:,:);
clear sample

% Reference
ref = load(fullfile(baseDir,fileRef));
medium.rfLR = resample(ref.rf1(:,:,:),p,q);
medium.rfHR = resample(ref.rf2(:,:,:),p,q);
medium.rfLR = medium.rfLR(1:length(medium.z),:,:);
medium.rfHR = medium.rfHR(1:length(medium.z),:,:);
clear ref

filterParams.freqC = freq;
[bz,xP,zP] = getMeasurements(medium,filterParams,blockParams);

idz = medium.z>0.5e-2;
bmodeFreq = db(hilbert(medium.rfH(idz,:,1)));
bmode(:,:,1) = bmodeFreq - max(bmodeFreq(:));
zBm = medium.z(idz); xBm = medium.x;

%% ATTENUATION ESTIMATION
% get Spectra
[Sp,Sd,x_ACS,z_ACS,f] = getSpectrum(medium.rfL,medium.x*100,medium.z*100, ...
    medium.fs,blockParamsAcs);
att_ref_map(1,1,:) = medium.alphaRcoeff/100*f.^medium.alphaRpower;
[SpRef,SdRef,~,~,~] = getSpectrum(medium.rfLR,medium.x*100,medium.z*100, ...
    medium.fs,blockParamsAcs);

% Setting up system
L = (z_ACS(2) - z_ACS(1))/(1 - blockParamsAcs.overlap)/2;   % (cm)
sld = log(Sp) - log(Sd);
sldRef = log(SpRef) - log(SdRef);
compensation = sldRef - 4*L*att_ref_map;
range = f>freqL/1e6 & f<freqH/1e6;
b = sld(:,:,range) - compensation(:,:,range);
ufr = f(range);
[m,n,p] = size(b);
A1 = kron( 4*L*ufr , speye(m*n) );
A2 = kron( ones(size(ufr)) , speye(m*n) );
mask = ones(m,n,p);

%% RSLD
tic
[Bn,Cn] = AlterOpti_ADMM(A1,A2,b(:),muRsld,muRsld,m,n,tol,mask(:));
toc
BR = (reshape(Bn,m,n));
CR = (reshape(Cn,m,n));

% Plot and cumulative map
attSam = (BR*filterParams.freqC)*100; % Np/cm to Np/m
% figure('Position',imPosition),
% imagesc(xP*100,zP*100,attSam)
% axis image
% c = colorbar;
% c.Label.String = "dB/cm";
% colormap turbo
% title('Attenuation from RSLD')
attSam = cumsum(attSam)./(1:length(zP))';
estBAdep = (bz./( 1-exp(-2*attSam.*zP) ) .*attSam.*zP - 1)*2 ;

%% Local maps with regularization
[m,n] = size(bz);
tol = 1e-3;
muLocal = 0.01;
dzP = zP(2)-zP(1);
izP = round(zP./dzP);
factorq = izP(1)./izP;
P = sparse(tril(ones(m-1)));
P = P./izP(2:end);
P = kron(speye(n),P);

estBAcum = estBAdep - estBAdep(1,:).*factorq;
estBAcum = estBAcum(2:end,:);
estBAinst = IRLS_TV(estBAcum(:),P,muLocal,m-1,n,tol,[],ones((m-1)*n,1));
estBAinst = reshape(estBAinst,m-1,n);

estAC = BR*db(exp(1));
%%
[X,Z] = meshgrid(xBm,zBm);
inc = (X-cx).^2./rx^2 + (Z-cz).^2./rz^2 < 0.9^2;
back = (X-cx).^2./rx^2 + (Z-cz).^2./rz^2 > 1.1^2;
[Xq,Zq] = meshgrid(xBm,zBm);
[X,Z] = meshgrid(x_ACS/100,z_ACS/100);
AttInterp = interp2(X,Z,estAC,Xq,Zq);
[X,Z] = meshgrid(xP,zP(2:end));
BaInterp = interp2(X,Z,estBAinst,Xq,Zq);

fprintf('AC inc: %.2f +/- %.2f\n', mean(AttInterp(inc),'omitnan'), ...
    std(AttInterp(inc), [] ,'omitnan'));
fprintf('AC back: %.2f +/- %.2f\n', mean(AttInterp(back),'omitnan'), ...
    std(AttInterp(back), [] ,'omitnan'));
fprintf('B/A inc: %.2f +/- %.2f\n', mean(BaInterp(inc),'omitnan'), ...
    std(BaInterp(inc), [] ,'omitnan'));
fprintf('B/A back: %.2f +/- %.2f\n', mean(BaInterp(back),'omitnan'), ...
    std(BaInterp(back), [] ,'omitnan'));

%%
[X,Z] = meshgrid(xBm,zBm);
incPlot = (X-cx).^2./rx^2 + (Z-cz).^2./rz^2 < 1^2;

roi = ones(size(bmode,[1 2]));
figure('Position',imPosition);
[~,hB,hColor] = imOverlayInterp(bmode(:,:,1),estBAinst,[-50 0],baRange,0.7,...
    xP*100,zP(2:end)*100,roi,xBm*100,zBm*100);
title('B/A')
xlabel('Lateral [cm]')
ylabel('Depth [cm]')
colormap pink
ylim(ylimBm)
hold on
contour(xBm*100,zBm*100,incPlot,1,'w--', 'LineWidth',2)
% contour(xBm*100,zBm*100,~isnan(BaInterp)&incPlot,1,'w--', 'LineWidth',2)
% contour(xBm*100,zBm*100,~isnan(BaInterp),1,'k--', 'LineWidth',2)
hold off

figure('Position',imPosition);
[~,hB,hColor] = imOverlayInterp(bmode(:,:,1),estAC,[-50 0],attRange,0.7,...
    x_ACS,z_ACS,roi,xBm*100,zBm*100);
title('\alpha_0')
xlabel('Lateral [cm]')
ylabel('Depth [cm]')
colormap turbo
hColor.Label.String = 'db/cm/MHz';
ylim(ylimBm)
hold on
contour(xBm*100,zBm*100,incPlot,1,'w--', 'LineWidth',2)
% contour(xBm*100,zBm*100,~isnan(BaInterp)&incPlot,1,'w--', 'LineWidth',2)
% contour(xBm*100,zBm*100,~isnan(BaInterp),1,'k--', 'LineWidth',2)
hold off

%% SWIFT
% First iteration
[~,Cn] = optimAdmmTvTikhonov(A1,A2,b(:),muBswift,muCswift,m,n,tol,mask(:));
bscMap = reshape(Cn*NptodB,m,n);

% Weight map
w = (1-reject)*(abs(bscMap)<ratioCutOff)+reject;
wExt = movmin(w,extension);

% Weight matrices and new system
W = repmat(wExt,[1 1 p]);
W = spdiags(W(:),0,m*n*p,m*n*p);
bw = W*b(:);        
A1w = W*A1;
A2w = W*A2;

% Second iteration
[Bn,cN] = optimAdmmWeightedTvTikhonov(A1w,A2w,bw,muBswift,muCswift,m,n,tol,mask(:),w);
BR = reshape(Bn,m,n);
CR = reshape(Cn,m,n);


% Plot and cumulative map
attSam = (BR*filterParams.freqC)*100; % Np/cm to Np/m
% figure('Position',imPosition),
% imagesc(xP*100,zP*100,attSam)
% axis image
% c = colorbar;
% c.Label.String = "dB/cm";
% colormap turbo
% title('Attenuation from SWIFT')
attSam = cumsum(attSam)./(1:length(zP))';
estBAdep = (bz./( 1-exp(-2*attSam.*zP) ) .*attSam.*zP - 1)*2 ;

%% Local maps with regularization
[m,n] = size(bz);
tol = 1e-3;
muLocal = 0.01;
dzP = zP(2)-zP(1);
izP = round(zP./dzP);
factorq = izP(1)./izP;
P = sparse(tril(ones(m-1)));
P = P./izP(2:end);
P = kron(speye(n),P);

estBAcum = estBAdep - estBAdep(1,:).*factorq;
estBAcum = estBAcum(2:end,:);
estBAinst = IRLS_TV(estBAcum(:),P,muLocal,m-1,n,tol,[],ones((m-1)*n,1));
estBAinst = reshape(estBAinst,m-1,n);

estAC = BR*db(exp(1));
%%
[X,Z] = meshgrid(xBm,zBm);
inc = (X-cx).^2./rx^2 + (Z-cz).^2./rz^2 < 0.9^2;
back = (X-cx).^2./rx^2 + (Z-cz).^2./rz^2 > 1.1^2;
[Xq,Zq] = meshgrid(xBm,zBm);
[X,Z] = meshgrid(x_ACS/100,z_ACS/100);
AttInterp = interp2(X,Z,estAC,Xq,Zq);
[X,Z] = meshgrid(xP,zP(2:end));
BaInterp = interp2(X,Z,estBAinst,Xq,Zq);

fprintf('AC inc: %.2f +/- %.2f\n', mean(AttInterp(inc),'omitnan'), ...
    std(AttInterp(inc), [] ,'omitnan'));
fprintf('AC back: %.2f +/- %.2f\n', mean(AttInterp(back),'omitnan'), ...
    std(AttInterp(back), [] ,'omitnan'));
fprintf('B/A inc: %.2f +/- %.2f\n', mean(BaInterp(inc),'omitnan'), ...
    std(BaInterp(inc), [] ,'omitnan'));
fprintf('B/A back: %.2f +/- %.2f\n', mean(BaInterp(back),'omitnan'), ...
    std(BaInterp(back), [] ,'omitnan'));

%%
[X,Z] = meshgrid(xBm,zBm);
incPlot = (X-cx).^2./rx^2 + (Z-cz).^2./rz^2 < 1^2;

roi = ones(size(bmode,[1 2]));
figure('Position',imPosition);
[~,hB,hColor] = imOverlayInterp(bmode(:,:,1),estBAinst,[-50 0],baRange,0.7,...
    xP*100,zP(2:end)*100,roi,xBm*100,zBm*100);
title('B/A')
xlabel('Lateral [cm]')
ylabel('Depth [cm]')
colormap pink
ylim(ylimBm)
hold on
contour(xBm*100,zBm*100,incPlot,1,'w--', 'LineWidth',2)
% contour(xBm*100,zBm*100,~isnan(BaInterp)&incPlot,1,'w--', 'LineWidth',2)
% contour(xBm*100,zBm*100,~isnan(BaInterp),1,'k--', 'LineWidth',2)
hold off

figure('Position',imPosition);
[~,hB,hColor] = imOverlayInterp(bmode(:,:,1),estAC,[-50 0],attRange,0.7,...
    x_ACS,z_ACS,roi,xBm*100,zBm*100);
title('\alpha_0')
xlabel('Lateral [cm]')
ylabel('Depth [cm]')
colormap turbo
hColor.Label.String = 'db/cm/MHz';
ylim(ylimBm)
hold on
contour(xBm*100,zBm*100,incPlot,1,'w--', 'LineWidth',2)
% contour(xBm*100,zBm*100,~isnan(BaInterp)&incPlot,1,'w--', 'LineWidth',2)
% contour(xBm*100,zBm*100,~isnan(BaInterp),1,'k--', 'LineWidth',2)
hold off

%%
save_all_figures_to_directory(resultsDir,'heteroFig','svg')
close all
