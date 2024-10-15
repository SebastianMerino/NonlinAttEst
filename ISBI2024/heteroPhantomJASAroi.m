% New methods with frequency compounding. Requires three transmissions
% Used for showing results

setup;
baseDir = 'C:\Users\sebas\Documents\Data\Nonlinearity\01-Oct-24\bf';
resultsDir = 'C:\Users\sebas\Documents\Data\Nonlinearity\resultsJASA\phantomTest';
[~,~,~] = mkdir(resultsDir);

% Auxiliar variables
NptodB = 20*log10(exp(1));

rz = 1.5e-2; rx = 1.3e-2;
cz = 4.2e-2; cx = 1.5e-2;

baRange = [4 12];
attRange = [0.05,0.2];
ylimBm = [1.5 5.8];
imPosition = [100 200 250 250];
alphaInit = 0.076;
baInit = 5.4;

c0 = 1470;
%% Preparing data
freqVec = [5,6,7]; % FRECUENCIES FOR FILTERING

% Known variables
medium.v = 5; % v= 5.438961611

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
blockParams.zlim = [2.5; 4.8]/100;
blockParams.xlim = [-1.5; 1.5]/100;
blockParams.downFactor = 20;


for iRoi = 1:5
    fprintf('\n\nROI %d\n',iRoi);
    switch iRoi
        case 1
            blockParams.zlim = [2.5; 4.8]/100;
            blockParams.xlim = [-1.5; 1.5]/100;
        case 2
            blockParams.zlim = [2.5; 5]/100;
            blockParams.xlim = [-1.5; 1.5]/100;
        case 3
            blockParams.zlim = [2.5; 5.5]/100;
            blockParams.xlim = [-1.5; 1.5]/100;
        case 4
            blockParams.zlim = [2; 4.8]/100;
            blockParams.xlim = [-1.5; 1.5]/100;
        case 5
            blockParams.zlim = [2.5; 4.8]/100;
            blockParams.xlim = [-1.8; 2]/100;

    end

%% Measurements, IUS version
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

% Measurements
[bz,zbz,xB,zB] = getMeasurementsIUS(medium,filterParams,blockParams);
X = permute(zbz,[3 1 2]);
Y = permute(bz,[3 1 2]);
X = X(:); Y = Y(:);

idz = medium.z>0.5e-2;
bmodeFreq = db(hilbert(medium.rfH(idz,:,1)));
bmode(:,:,1) = bmodeFreq - max(bmodeFreq(:));
zBm = medium.z(idz); xBm = medium.x;

%% Gauss-Newton with LM
[m,n,p] = size(bz);

tol = 1e-3;
maxIte = 200;
muAlpha = 0.1; muBeta = 0.1;
beta0 = 1+baInit/2;
alpha0 = alphaInit*freq^2 *100/NptodB; % Np/m

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

estAClm = reshape(alphaArr/freq^2 /100*NptodB,[m,n]);
estBAlm = reshape(2*(betaArr-1),[m,n]);

%%
[X,Z] = meshgrid(xBm,zBm);
inc = (X-cx).^2./rx^2 + (Z-cz).^2./rz^2 < 0.9^2;
back = (X-cx).^2./rx^2 + (Z-cz).^2./rz^2 > 1.1^2;
[X,Z] = meshgrid(xB,zB);
[Xq,Zq] = meshgrid(xBm,zBm);
AttInterp = interp2(X,Z,estAClm,Xq,Zq);
BaInterp = interp2(X,Z,estBAlm,Xq,Zq);

fprintf('AC inc: %.2f +/- %.2f\n', mean(AttInterp(inc),'omitnan'), ...
    std(AttInterp(inc), [] ,'omitnan'));
fprintf('AC back: %.2f +/- %.2f\n', mean(AttInterp(back),'omitnan'), ...
    std(AttInterp(back), [] ,'omitnan'));
fprintf('B/A inc: %.2f +/- %.2f\n', mean(BaInterp(inc),'omitnan'), ...
    std(BaInterp(:), [] ,'omitnan'));
fprintf('B/A back: %.2f +/- %.2f\n', mean(BaInterp(back),'omitnan'), ...
    std(BaInterp(:), [] ,'omitnan'));

%%
[X,Z] = meshgrid(xBm,zBm);
incPlot = (X-cx).^2./rx^2 + (Z-cz).^2./rz^2 < 1^2;

roi = ones(size(bmode,[1 2]));
figure('Position',imPosition);
[~,hB,hColor] = imOverlayInterp(bmode(:,:,1),estBAlm,[-50 0],baRange,0.7,...
    xB*100,zB*100,roi,xBm*100,zBm*100);
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
[~,hB,hColor] = imOverlayInterp(bmode(:,:,1),estAClm,[-50 0],attRange,0.7,...
    xB*100,zB*100,roi,xBm*100,zBm*100);
title('\alpha_0')
xlabel('Lateral [cm]')
ylabel('Depth [cm]')
colormap turbo
hColor.Label.String = 'db/cm/MHz^2';
ylim(ylimBm)
hold on
contour(xBm*100,zBm*100,incPlot,1,'w--', 'LineWidth',2)
% contour(xBm*100,zBm*100,~isnan(BaInterp)&incPlot,1,'w--', 'LineWidth',2)
% contour(xBm*100,zBm*100,~isnan(BaInterp),1,'k--', 'LineWidth',2)
hold off


%% Measurements, new version
medium.fs = 30e6; % new sampling frequency
medium.maxDepth = 6e-2;
dz = (1/medium.fs)*c0/2;
medium.z = (0:dz:medium.maxDepth)';

bzf = [];
for iFreq = 1:length(freqVec)
    %%
    freq = freqVec(iFreq);
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

    filterParams.freqC = freqVec(iFreq);
    [bz,xP,zP] = getMeasurements(medium,filterParams,blockParams);
    bzf(:,:,iFreq) = bz;

    idz = medium.z>0.5e-2;
    bmodeFreq = db(hilbert(medium.rfH(idz,:,1)));
    bmode(:,:,iFreq) = bmodeFreq - max(bmodeFreq(:));
    zBm = medium.z(idz); xBm = medium.x;

end

%% ADMM
% Initialization
[Xmesh,Zmesh] = meshgrid(xP,zP);
alphaL = ones(size(Xmesh))*alphaInit;
betaL = ones(size(Xmesh))*(1+baInit/2);
u0 = [alphaL(:)*100/NptodB;betaL(:)];

% Optimizes F(u) + R(v)
% Hyperparameters
[m,n,p] = size(bzf);
tol = 1e-4;
muAlpha = 0.01; muBeta = 0.001; % testing 03-10
rho = 0.1;
maxIte = 300;

% Initialization
u = u0;
v = u;
w = zeros(size(u));

% Auxiliar
dzP = zP(2)-zP(1);
izP = round(zP./dzP);
DP = [-izP izP];
DP = spdiags(DP, [-1 0],m,m);
DP(1,:) = DP(2,:);
DP = kron(speye(n),DP);
mask = ones(n*m,1);
I = speye(2*m*n);
Ii = I(:,1:m*n);
Id = I(:,1+m*n:end);

% Objective functions
Fid = []; Reg = []; Dual = [];
Fid(1) = 1/2*norm( modelFreq(u,zP,freqVec) - bzf(:) )^2;
Reg(1) = muAlpha*TVcalc_isotropic(DP*u(1:m*n),m,n,mask) + ...
    muBeta*TVcalc_isotropic(DP*u(m*n+1:end),m,n,mask);
Dual(1) = 0;
ite  = 0;
error = 1;
% figure; tiledlayout(1,2)
% t1 = nexttile;         t2 = nexttile;
tic
while ite < maxIte % abs(error) > tol &&
    ite = ite + 1;

    % Fidelity step
    u = optimNonLinearGNFreq(zP,freqVec,bzf(:), speye(2*m*n),v-w, rho, tol,u);

    % Regularization step
    v = optimAdmmTvDy(Ii,Id,u+w, muAlpha/rho,muBeta/rho ,m,n,tol,mask,DP);

    % Dual update
    w = w + u - v;

    % Loss
    Fid(ite+1) = 1/2*norm( modelFreq(u,zP,freqVec) - bzf(:) )^2;
    Reg(ite+1) = muAlpha*TVcalc_isotropic(DP*v(1:m*n),m,n,mask) + ...
        muBeta*TVcalc_isotropic(DP*v(m*n+1:end),m,n,mask);
    Dual(ite+1) = norm(u-v);
    % fprintf("Ite: %01i, Fid: %.3f, Reg: %.3f, Dual: %.3f\n",...
    %     ite,Fid(ite+1),Reg(ite+1),Dual(ite+1))
    error = Fid(ite+1) + Reg(ite+1) - Fid(ite) - Reg(ite);

end
toc
%%

alphaArr = reshape(DP*u(1:m*n),[m,n]);
betaArr = reshape(DP*u(1+m*n:end),[m,n]);

estACtv = alphaArr*NptodB/100;
estBAtv = 2*(betaArr-1);

%%
[X,Z] = meshgrid(xBm,zBm);
inc = (X-cx).^2./rx^2 + (Z-cz).^2./rz^2 < 0.9^2;
back = (X-cx).^2./rx^2 + (Z-cz).^2./rz^2 > 1.1^2;
[X,Z] = meshgrid(xP,zP);
[Xq,Zq] = meshgrid(xBm,zBm);
AttInterp = interp2(X,Z,estACtv,Xq,Zq);
BaInterp = interp2(X,Z,estBAtv,Xq,Zq);

fprintf('AC inc: %.2f +/- %.2f\n', mean(AttInterp(inc),'omitnan'), ...
    std(AttInterp(inc), [] ,'omitnan'));
fprintf('AC back: %.2f +/- %.2f\n', mean(AttInterp(back),'omitnan'), ...
    std(AttInterp(back), [] ,'omitnan'));
fprintf('B/A inc: %.2f +/- %.2f\n', mean(BaInterp(inc),'omitnan'), ...
    std(BaInterp(:), [] ,'omitnan'));
fprintf('B/A back: %.2f +/- %.2f\n', mean(BaInterp(back),'omitnan'), ...
    std(BaInterp(:), [] ,'omitnan'));
%%
[X,Z] = meshgrid(xBm,zBm);
incPlot = (X-cx).^2./rx^2 + (Z-cz).^2./rz^2 < 1^2;

roi = ones(size(bmode,[1 2]));
figure('Position',imPosition);
[~,hB,hColor] = imOverlayInterp(bmode(:,:,2),estBAtv,[-50 0],baRange,0.7,...
    xP*100,zP*100,roi,xBm*100,zBm*100);
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
[~,hB,hColor] = imOverlayInterp(bmode(:,:,2),estACtv,[-50 0],attRange,0.7,...
    xP*100,zP*100,roi,xBm*100,zBm*100);
title('\alpha_0')
xlabel('Lateral [cm]')
ylabel('Depth [cm]')
colormap turbo
hColor.Label.String = '[db/cm/MHz^2]';
ylim(ylimBm)
hold on
contour(xBm*100,zBm*100,incPlot,1,'w--', 'LineWidth',2)
% contour(xBm*100,zBm*100,~isnan(BaInterp)&incPlot,1,'w--', 'LineWidth',2)
% contour(xBm*100,zBm*100,~isnan(BaInterp),1,'k--', 'LineWidth',2)
hold off
%%
save_all_figures_to_directory(resultsDir,char("heteroRoi"+iRoi+"Fig"),'svg')
close all


end
% ---------------------------------------------------------------------- %
% ---------------------------------------------------------------------- %
%%
function [bz,xP,zP] = getMeasurements(medium,filterParams,blockParams)
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

% Block averaging
[meanP,xP,zP] = getMeanBlock(cat(3,PLfull,PHfull,PLRfull,PHRfull),x,z,...
    blockParams);
PL = meanP(:,:,1);
PH = meanP(:,:,2);
PLR = meanP(:,:,3);
PHR = meanP(:,:,4);

% Getting measurements
attRef = medium.alphaRcoeff*(freqC/1e6)^medium.alphaRpower;
% attRef = medium.alphaR*(freqC/1e6)^2;
bz = medium.betaR*sqrt( abs(v*PL-PH)./abs(v*PLR-PHR) .*PLR./PL ).*...
    ( 1-exp(-2*attRef*zP) )/attRef./zP;

end
