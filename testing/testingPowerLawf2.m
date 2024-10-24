% Genera imágenes de data simulada uniforme
% Referencia con f2, sample con f1p2, modelo con f2
setup;
baseDir = 'P:\smerino\NonlinAttESt\data\ref_att1p2';
refDir = 'P:\smerino\NonlinAttESt\data\uniformSimulation';
resultsDir = fullfile(baseDir,'results_BA10_att0p60f1p2');
[~,~,~] = mkdir(resultsDir);

% Cambiar para seleccionar muestra
for iSim = 1:8   % del 1 al 8

%% Auxiliar variables
NptodB = 20*log10(exp(1));
radiusDisk = (9)*1e-3;
centerDepth = 22.5e-3;

baRange = [5 13];
attRange = [0.05,0.15];
imPosition = [100 200 250 300];

alphaRefVec = 0.08:0.02:0.22;

%% Preparing data
baInit = 6;

% Known variables
medium.v = 5;
medium.betaR = 1 + 6/2;

% Filtering parameters
filterParams.freqC = 5;
filterParams.freqTol = 0.5;
filterParams.nCycles = 10; % Number of cycles of the initial filter

% Subsampling parameters
wl = 1540/5e6; % Mean central frequency
blockParams.blockSize = [25 25]*wl;
blockParams.overlap = 0.8;
blockParams.zlim = [1; 5.5]/100;
blockParams.xlim = [-2.5; 2.5]/100;
blockParams.downFactor = 20;

freqVec = [4,5,6]; % FRECUENCIES FOR FILTERING

%% For loop
alphaRef = alphaRefVec(iSim);
alphaInit = alphaRef;
medium.alphaR = alphaRef/NptodB*100; % alpha0 in dB/100/MHz2

%% Measurements, IUS version
freq = 5;
fileSam = "RFfn2_PWNE"+freq+"MHz_refBA10_att0p60f1p2_nc10_400kPa";
fileRef = "RFfn2_PWNE"+freq+"MHz_refBA6_att0p"+...
    num2str(round(100*alphaRef), '%02d')+"f2p0_400kPa";

% Sample
sample = load(fullfile(baseDir,fileSam));
medium.z = sample.z';
medium.x = sample.x;
medium.fs = sample.fs;
medium.rfL = sample.rf1(:,:,:);
medium.rfH = sample.rf2(:,:,:);
clear sample

% Reference
ref = load(fullfile(refDir,fileRef));
medium.rfLR = ref.rf1(:,:,:);
medium.rfHR = ref.rf2(:,:,:);
clear ref

% Measurements
[bz,zbz,xB,zB] = getMeasurementsIUS(medium,filterParams,blockParams);
X = permute(zbz,[3 1 2]);
Y = permute(bz,[3 1 2]);
X = X(:); Y = Y(:);

%% Gauss-Newton with LM
[m,n,p] = size(bz);

tol = 1e-3;
maxIte = 200;
% muAlpha = 1000; muBeta = 0;
muAlpha = 1; muBeta = 1;
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

baMean = mean(estBAlm(:),'omitnan');
baStd = std(estBAlm(:),[],'omitnan');
acMean = mean(estAClm(:),'omitnan');
acStd = std(estAClm(:),[],'omitnan');

figure('Position',imPosition);
imagesc(xB*1e2,zB*1e2,estBAlm); colorbar;
clim(baRange);
title("B/A = "+num2str(baMean,2)+"\pm"+num2str(baStd,2));
axis image
colormap pink; colorbar;
xlabel('Lateral [cm]');
ylabel('Depth [cm]');

figure('Position',imPosition);
imagesc(xB*1e2,zB*1e2,estAClm); colorbar;
clim(attRange);
title("\alpha_0 = "+num2str(acMean,2)+"\pm"+num2str(acStd,2));
axis image
colormap turbo; colorbar;
xlabel('Lateral [cm]');
ylabel('Depth [cm]');
pause(0.1)

metricsIUS(iSim) = getMetrics(estAClm,estBAlm,'IUS',alphaRef);


%% Measurements, new version
bzf = [];
for iFreq = 1:length(freqVec)
    freq = freqVec(iFreq);
    fileSam = "RFfn2_PWNE"+freq+"MHz_refBA10_att0p60f1p2_nc10_400kPa";
    fileRef = "RFfn2_PWNE"+freq+"MHz_refBA6_att0p"+...
        num2str(round(100*alphaRef), '%02d')+"f2p0_400kPa";

    % Sample
    sample = load(fullfile(baseDir,fileSam));
    medium.z = sample.z';
    medium.x = sample.x;
    medium.fs = sample.fs;
    medium.rfL = sample.rf1(:,:,:);
    medium.rfH = sample.rf2(:,:,:);
    clear sample

    % Reference
    ref = load(fullfile(refDir,fileRef));
    medium.rfLR = ref.rf1(:,:,:);
    medium.rfHR = ref.rf2(:,:,:);
    clear ref

    filterParams.freqC = freqVec(iFreq);
    [bz,xP,zP] = getMeasurements(medium,filterParams,blockParams);
    bzf(:,:,iFreq) = bz;
end

% Initialization
[Xmesh,Zmesh] = meshgrid(xP,zP);
alphaL = ones(size(Xmesh))*alphaInit;
betaL = ones(size(Xmesh))*(1+baInit/2);
u0 = [alphaL(:)*100/NptodB;betaL(:)];


%% ADMM
% Optimizes F(u) + R(v)
% Hyperparameters
[m,n,p] = size(bzf);
tol = 1e-3;
muAlpha = 1; muBeta = 0.1;
rho = 1;
maxIte = 200;

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
tic
while abs(error) > tol && ite < maxIte
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
    fprintf("Ite: %01i, Fid: %.3f, Reg: %.3f, Dual: %.3f\n",...
        ite,Fid(ite+1),Reg(ite+1),Dual(ite+1))
    error = Fid(ite+1) + Reg(ite+1) - Fid(ite) - Reg(ite);
end
toc
%%
alphaArr = reshape(DP*u(1:m*n),[m,n]);
betaArr = reshape(DP*u(1+m*n:end),[m,n]);

estACtv = alphaArr*NptodB/100;
estBAtv = 2*(betaArr-1);

baMean = mean(estBAtv(:),'omitnan');
baStd = std(estBAtv(:),[],'omitnan');
acMean = mean(estACtv(:),'omitnan');
acStd = std(estACtv(:),[],'omitnan');

figure('Position',imPosition);
imagesc(xP*1e3,zP*1e3,estBAtv); colorbar;
clim(baRange);
title("B/A = "+num2str(baMean,2)+"\pm"+num2str(baStd,2));
axis image
colormap pink; colorbar;
xlabel('Lateral [cm]');
ylabel('Depth [cm]');

figure('Position',imPosition);
im = imagesc(xP*1e3,zP*1e3,estACtv); colorbar;
clim(attRange);
title("\alpha_0 = "+num2str(acMean,2)+"\pm"+num2str(acStd,2));
axis image
colormap turbo; colorbar;
xlabel('Lateral [cm]');
ylabel('Depth [cm]');
pause(0.1)

metricsADMM(iSim) = getMetrics(estACtv,estBAtv,'ADMM',alphaRef);

%%
save_all_figures_to_directory(resultsDir,...
    char("sim"+iSim+"fig"))
close all
end
%%
disp(struct2table(metricsIUS(iSim)))
disp(struct2table(metricsADMM(iSim)));
T = [struct2table(metricsIUS);struct2table(metricsADMM)];
writetable(T,fullfile(resultsDir,'homoRef.xlsx'))

% ---------------------------------------------------------------------- %
% ---------------------------------------------------------------------- %
%%
function metrics = getMetrics(ac,ba,method,alphaRef)
metrics.acMean = mean(ac(:),'omitnan');
metrics.acStd = std(ac(:),[],'omitnan');
metrics.baMean = mean(ba(:),'omitnan');
metrics.baStd = std(ba(:),[],'omitnan');
metrics.method = method;
metrics.alphaRef = alphaRef;
end