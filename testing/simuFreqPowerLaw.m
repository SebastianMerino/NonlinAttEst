setup;
baseDir = 'P:\smerino\simulation_acs\rf_data\24_10_28_ba2d\fn2';
refDir = 'P:\smerino\simulation_acs\rf_data\24_10_28_ba2d\fn2';
resultsDir = 'P:\smerino\simulation_acs\rf_data\24_10_28_ba2d\results';
[~,~,~] = mkdir(resultsDir);

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
alphaInit = 0.4;
gammaInit = 1.2;

% Known variables
medium.v = 5;
medium.betaR = 1 + 6/2;
medium.alphaR = 0.4/NptodB*100; % alpha0 in dB/100/MHz2
gammaRef = 1.2;

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

freqVec = 4:0.5:8; % FRECUENCIES FOR FILTERING

%% Measurements, new version
bzf = [];
for iFreq = 1:length(freqVec)
    freq = freqVec(iFreq);
    disp(freq)
    fileRef = "RFfn2_PWNE_"+ strrep(num2str(freq),'.','p') + ...
        "MHz_homo_BA6_att0p4f1p2_nc10_400kPa";
    fileSam = "RFfn2_PWNE_"+ strrep(num2str(freq),'.','p') + ...
        "MHz_homo_BA10_att0p6f1p2_nc10_400kPa";

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
    [bz,xP,zP] = getMeasurements_v2(medium,filterParams,blockParams,gammaRef);
    bzf(:,:,iFreq) = bz;
end

% Initialization
[Xmesh,Zmesh] = meshgrid(xP,zP);
alphaL = ones(size(Xmesh))*alphaInit;
betaL = ones(size(Xmesh))*(1+baInit/2);
u0 = [alphaL(:)*100/NptodB;betaL(:)];

%% Gauss-Newton with LM
[m,n,p] = size(bzf);
tol = 1e-4;
maxIte = 20;
% betaIni = 1+baInit/2;
% alphaIni = alphaInit/NptodB*100;
% gammaIni = gammaInit;
betaIni = 1+10/2;
alphaIni = 0.6/NptodB*100;
gammaIni = 1.2;

muAlpha = 0; muBeta = 0; muGamma = 1000;

u = [alphaIni*ones(n*m,1);gammaIni*ones(n*m,1);betaIni*ones(n*m,1)];
regMatrix = blkdiag(muAlpha*speye(n*m),muGamma*speye(m*n),...
    muBeta*speye(n*m));

loss = [];
loss(1) = 0.5*norm(modelFreqPL(u,zP,freqVec) - bzf(:))^2;

ite = 1;
tic 
while true
    jcb = jacobianFreqPL(u,zP,freqVec);
    res = modelFreqPL(u,zP,freqVec) - bzf(:);
    [step,~] = pcg(jcb'*jcb + regMatrix,jcb'*-res, tol, 200);
    u = u + step;

    loss(ite+1) = 0.5*norm(modelFreqPL(u,zP,freqVec) - bzf(:))^2;
    if abs(loss(ite+1) - loss(ite))<tol || ite == maxIte, break; end
    ite = ite + 1;
end
toc
%%
alpha0Map = reshape(u(1:m*n),m,n);
gammaMap = reshape(u(m*n+1:2*m*n),m,n);
betaMap = reshape(u(2*m*n+1:end),m,n);

estAClm = reshape(alpha0Map*NptodB/100,[m,n]);
estBAlm = reshape(2*(betaMap-1),[m,n]);

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
clim(attRange);
title('\alpha_0 in \alpha(f) = \alpha_0 \times f^2 dB/100')
axis image
colormap turbo; colorbar;
xlabel('Lateral distance (mm)');
ylabel('Depth (mm)');

figure; imagesc(xP*1e3,zP*1e3,gammaMap); colorbar;
clim([1 2.1]);
title('\gamma')
axis image
colormap turbo; colorbar;
xlabel('Lateral distance (mm)');
ylabel('Depth (mm)');


figure; imagesc(xP*1e3,zP*1e3,estBAlm); colorbar;
clim(baRange);
title('B/A');
title(['B/A = ',num2str(median(estBAlm(:),'omitnan'),'%.1f')]);
axis image
colormap pink; colorbar;
xlabel('Lateral distance (mm)');
ylabel('Depth (mm)');
%title("\bf BoA map");

pause(0.5)
