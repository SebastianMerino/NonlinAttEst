% New methods with frequency compounding. Requires three transmissions
% Used for showing results

setup;
% baseDir = 'C:\Users\sebas\Documents\Data\Nonlinearity\homoAttMismatch\f2';
% baseDir = 'C:\Users\sebas\Documents\Data\Nonlinearity\rfdata';
baseDir = 'C:\Users\sebas\Documents\Data\Nonlinearity\attInc\fn2';

% Auxiliar variables
NptodB = 20*log10(exp(1));
radiusDisk = (9)*1e-3;
centerDepth = 22.5e-3;

baRange = [5 13];
attRange = [0.09,0.17];

alphaInc = 0.16;

alpha0Init = 0;
alpha1Init = 0.5;
betaInit = 6;

%% Preparing data
% Known variables
medium.v = 5;
medium.betaR = 1 + 12/2;             
medium.alphaR = 0.1/NptodB*100; % alpha0 in dB/100/MHz2

% Filtering parameters
filterParams.freqC = 5;
filterParams.freqTol = 0.5;
filterParams.nCycles = 10; % Number of cycles of the initial filter

% Subsampling parameters
wl = 1540/5e6; % Mean central frequency
blockParams.blockSize = [20 20]*wl; 
blockParams.overlap = 0.8;
blockParams.zlim = [0.5; 5.5]/100;
blockParams.xlim = [-2.5; 2.5]/100;

freqVec = [4,5,6]; % FRECUENCIES FOR FILTERING

%% Getting B/A
bzf = [];
for iFreq = 1:length(freqVec)
    freq = freqVec(iFreq);
    alphaStr = num2str(alphaInc*100,"%02d");
    fileSam = "RFfn2_PWNE"+freq+"MHz_samincBA6inc12_att0p1f2inc0p"+alphaStr+ ...
        "_nc10_400kPa";
    fileRef = "RFfn2_PWNE"+freq+"MHz_samBA12_att0p1f2inc0p10_nc10_400kPa";

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
    [bz,xP,zP] = getMeasurements(medium,filterParams,blockParams);
    bzf(:,:,iFreq) = bz;

    bmode = db(hilbert(medium.rfH(300:end-600,:,1)));
    bmode = bmode - max(bmode (:));
    figure,
    imagesc(medium.x*100,medium.z(300:end-600)*100,bmode, [-40 0])
    colormap gray
    axis image
    colorbar
    title('B-mode')
end

% Measurements
figure('Units','centimeters', 'Position',[5 5 20 6]),
tiledlayout(1,3)
nexttile,
imagesc(xP*100,zP*100,bzf(:,:,1), [2 10])
title("Measurements b(z,f)")
xlabel('Lateral [cm]')
ylabel('Depth [cm]')
axis image
colorbar
colormap parula

nexttile,
imagesc(xP*100,zP*100,bzf(:,:,2), [2 10])
title("Measurements b(z,f)")
xlabel('Lateral [cm]')
ylabel('Depth [cm]')
axis image
colorbar
colormap parula

nexttile,
imagesc(xP*100,zP*100,bzf(:,:,3), [2 10])
title("Measurements b(z,f)")
xlabel('Lateral [cm]')
ylabel('Depth [cm]')
axis image
colorbar
colormap parula

%% Initialization
[Xmesh,Zmesh] = meshgrid(xP,zP);
alpha0L = ones(size(Xmesh))*alpha0Init;
gamma = ones(size(Xmesh))*2;
betaL = ones(size(Xmesh))*(1+betaInit/2);
u0 = [alpha0L(:)*100/NptodB;gamma(:);betaL(:)];


%% Gauss-Newton with LM
[m,n,p] = size(bzf);
tol = 1e-4;
maxIte = 400;
betaIni = 1+(10.5)/2;
alphaIni = 0.08/NptodB*100;
gammaIni = 2;
muAlpha = 0; muBeta = 0; muGamma = 0;

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


%% Inversion, Coila's implementation
zBlock = zP; xBLock = xP; 
muLocal = 0.01;
dzBlock = zBlock(2)-zBlock(1);
izBlock = round(zBlock./dzBlock);

factorq = izBlock(1)./izBlock ;

estBAcum = estBAlm - estBAlm(1,:).*factorq; 
estBAcum = estBAcum(2:end,:);

estACcum = estAClm - estAClm(1,:).*factorq; 
estACcum = estACcum(2:end,:);


P = sparse(tril(ones(m-1)));
P = P./izBlock(2:end);
P = kron(speye(n),P);
estBAinst = IRLS_TV(estBAcum(:),P,muLocal,m-1,n,tol,[],ones((m-1)*n,1));
estBAinst = reshape(estBAinst,m-1,n);

estACinst = IRLS_TV(estACcum(:),P,muLocal/100,m-1,n,tol,[],ones((m-1)*n,1));
estACinst = reshape(estACinst,m-1,n);

figure; imagesc(xP*1e3,zP*1e3,estACinst); colorbar;
clim(attRange);
title('\alpha_0 in \alpha(f) = \alpha_0 \times f^2 dB/cm')
axis image
colormap turbo; colorbar;
xlabel('Lateral distance (mm)');
ylabel('Depth (mm)');
hold on
rectangle('Position',[0-radiusDisk,centerDepth-radiusDisk,...
2*radiusDisk,2*radiusDisk]*1000, 'Curvature',1,...
'EdgeColor','b', 'LineStyle','--', 'LineWidth',2)
hold off


figure; imagesc(xP*1e3,zP*1e3,estBAinst); colorbar;
clim(baRange);
title('B/A');
axis image
colormap pink; colorbar;
xlabel('Lateral distance (mm)');
ylabel('Depth (mm)');
hold on
rectangle('Position',[0-radiusDisk,centerDepth-radiusDisk,...
2*radiusDisk,2*radiusDisk]*1000, 'Curvature',1,...
'EdgeColor','b', 'LineStyle','--', 'LineWidth',2)
hold off
pause(0.1)



[Xmesh,Zmesh] = meshgrid(xP,zP(2:end));
inc = Xmesh.^2 + (Zmesh-centerDepth).^2 < (radiusDisk-3e-3).^2;
back = Xmesh.^2 + (Zmesh-centerDepth).^2 > (radiusDisk+3e-3).^2;
fprintf('AC inc: %.2f +/- %.2f\n', mean(estACinst(inc),'omitnan'), ...
std(estBAinst(inc), [] ,'omitnan'));
fprintf('AC back: %.2f +/- %.2f\n', mean(estACinst(back),'omitnan'), ...
std(estBAinst(back), [] ,'omitnan'));

fprintf('B/A inc: %.2f +/- %.2f\n', mean(estBAinst(inc),'omitnan'), ...
std(estBAinst(inc), [] ,'omitnan'));
fprintf('B/A back: %.2f +/- %.2f\n', mean(estBAinst(back),'omitnan'), ...
std(estBAinst(back), [] ,'omitnan'));