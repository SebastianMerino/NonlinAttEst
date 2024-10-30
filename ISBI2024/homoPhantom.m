% New methods with frequency compounding. Requires three transmissions
% Used for showing results

setup;
baseDir = 'C:\Users\sebas\Documents\Data\Nonlinearity\01-Oct-24\bf';

% Auxiliar variables
NptodB = 20*log10(exp(1));
radiusDisk = (9)*1e-3;
centerDepth = 22.5e-3;

baRange = [4 12];
attRange = [0.05,0.2];


alphaInit = 0.09;
baInit = 10.5;

c0 = 1470;
%% Preparing data
freqVec = [5,6,7]; % FRECUENCIES FOR FILTERING

% Known variables
medium.v = 5;

% % Reference adriana
medium.betaR = 1 + 5.4/2;
medium.alphaRcoeff = 0.3504/NptodB*100; %alpha0 in dB/m/MHz^2
medium.alphaRpower = 1.1483;

% MY REFERENCE

medium.fs = 30e6; % new sampling frequency
medium.maxDepth = 6e-2;
dz = (1/medium.fs)*c0/2;
medium.z = (0:dz:medium.maxDepth)';

% Filtering parameters
filterParams.freqC = mean(freqVec);
filterParams.freqTol = 0.5;
filterParams.nCycles = 10; % Number of cycles of the initial filter

% Subsampling parameters
wl = c0/filterParams.freqC/1e6; % Mean central frequency
blockParams.blockSize = [30 30]*wl;
blockParams.overlap = 0.8;
blockParams.zlim = [2.5; 5.5]/100;
blockParams.xlim = [-1.5; 1.5]/100;


%% Getting B/A
bzf = [];
for iFreq = 1:length(freqVec)
    %%
    freq = freqVec(iFreq);
    fileSam = "RFfn2_L11-5v_8v40v_"+freq+"MHz_uniform_cornstachx3_1";
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

    idz = medium.z>blockParams.zlim(1) & medium.z<blockParams.zlim(2);
    bmode = db(hilbert(medium.rfH(idz,:,1)));
    bmode = bmode - max(bmode (:));
    figure,
    imagesc(medium.x*100,medium.z(idz)*100,bmode, [-40 0])
    colormap gray
    axis image
    colorbar
    title('B-mode')
end

%% Measurements
% bLim = [2,10];
bLim = [2,4.5];

figure('Units','centimeters', 'Position',[5 5 20 6]),
tiledlayout(1,3)
nexttile,
imagesc(xP*100,zP*100,bzf(:,:,1), bLim)
title("Measurements b(z,f)")
xlabel('Lateral [cm]')
ylabel('Depth [cm]')
axis image
colorbar
colormap parula

nexttile,
imagesc(xP*100,zP*100,bzf(:,:,2), bLim)
title("Measurements b(z,f)")
xlabel('Lateral [cm]')
ylabel('Depth [cm]')
axis image
colorbar
colormap parula

nexttile,
imagesc(xP*100,zP*100,bzf(:,:,3), bLim)
title("Measurements b(z,f)")
xlabel('Lateral [cm]')
ylabel('Depth [cm]')
axis image
colorbar
colormap parula

% nexttile,
% imagesc(xP*100,zP*100,bzf(:,:,4), [2 7])
% title("Measurements b(z,f)")
% xlabel('Lateral [cm]')
% ylabel('Depth [cm]')
% axis image
% colorbar
% colormap parula

% %%
% save("testing\bzf_unix3_2_ref_2","bzf")
% out1 = load('testing\bzf_unix3_1_ref_3.mat');
% out2 = load('testing\bzf_unix3_2_ref_2.mat');
% bzf = (out1.bzf + out2.bzf)/2;
%% Initialization
[Xmesh,Zmesh] = meshgrid(xP,zP);
alphaL = ones(size(Xmesh))*alphaInit;
betaL = ones(size(Xmesh))*(1+baInit/2);
u0 = [alphaL(:)*100/NptodB;betaL(:)];

% inc = Xmesh.^2 + (Zmesh-centerDepth).^2 <= (radiusDisk).^2;
% [m,n,p] = size(bzf);
% alphaL(inc) = 0.16;
% izBlock = 1:length(zP);
% P = sparse(tril(ones(m)));
% P = P./izBlock';
% betaC = P*(betaL);
% alphaC =  P*(alphaL);
% u0 = [alphaC(:)*100/NptodB;betaC(:)];
%
% figure,
% tiledlayout(1,2)
% nexttile,
% imagesc(xP*1e3,zP*1e3,alphaC)
% clim(attRange);
% title('Init \alpha_0 in \alpha(f) = \alpha_0 \times f^2')
% axis image
% colormap turbo; colorbar;
% xlabel('Lateral distance (mm)');
% ylabel('Depth (mm)');
% t2 = nexttile;
% imagesc(xP*1e3,zP*1e3,(betaC-1)*2)
% clim(baRange);
% title('Init B/A');
% axis image
% colormap(t2,pink); colorbar;
% xlabel('Lateral distance (mm)');
% ylabel('Depth (mm)');
% pause(0.5)

%% Gauss-Newton with LM
[m,n,p] = size(bzf);
tol = 1e-3;
maxIte = 100;
muAlpha = 0;
muBeta = 0;
regMatrix = blkdiag(muAlpha*speye(n*m),muBeta*speye(n*m));

u = u0;
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

    alphaArr = u(1:n*m);
    betaArr = u(n*m+1:end);

    estAClm = reshape(alphaArr*NptodB/100,[m,n]);
    estBAlm = reshape(2*(betaArr-1),[m,n]);

end
toc

alphaArr = u(1:n*m);
betaArr = u(n*m+1:end);

estAClm = reshape(alphaArr*NptodB/100,[m,n]);
estBAlm = reshape(2*(betaArr-1),[m,n]);

fprintf('AC: %.3f +/- %.3f\n', mean(estAClm(:), 'omitnan'), ...
    std(estAClm(:), [] ,'omitnan'));
fprintf('B/A : %.2f +/- %.2f\n', mean(estBAlm(:),'omitnan'), ...
    std(estBAlm(:), [] ,'omitnan'));

figure('Units','centimeters', 'Position',[5 5 10 5]),
plot(loss, 'LineWidth',2)
% xlim([2 length(loss)])
xlabel('Number of iterations')
ylabel('Loss')


figure; imagesc(xP*1e3,zP*1e3,estAClm); colorbar;
clim(attRange);
title('\alpha_0 in \alpha(f) = \alpha_0 \times f^2 dB/100')
axis image
colormap turbo; colorbar;
xlabel('Lateral distance (mm)');
ylabel('Depth (mm)');

figure; imagesc(xP*1e3,zP*1e3,estBAlm); colorbar;
clim(baRange);
title('B/A');
axis image
colormap pink; colorbar;
xlabel('Lateral distance (mm)');
ylabel('Depth (mm)');
pause(0.5)


%% ADMM
% Optimizes F(u) + R(v)
% Hyperparameters
[m,n,p] = size(bzf);
tol = 1e-4;
muAlpha = 0.1; muBeta = 0.01; % homogeneous
rho = 0.1;
maxIte = 400;

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
% f = figure; tiledlayout(1,2)
% t1 = nexttile;         t2 = nexttile;
tic
while ite < maxIte %abs(error) > tol && 
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
    
    if mod(ite,10) ==1 % && ite<50
        alphaArr = reshape(DP*u(1:m*n),[m,n]);
        betaArr = reshape(DP*u(1+m*n:end),[m,n]);
        estACtv = alphaArr*NptodB/100;
        estBAtv = 2*(betaArr-1);
    end
end
toc
%%
close all,
figure('Units','centimeters', 'Position',[5 5 10 15]),
tiledlayout(2,1)
nexttile,
plot(Fid+Reg, 'LineWidth',2)
ylabel('Objective')
xlim([3 length(Fid)])
nexttile,
plot(Dual, 'LineWidth',2)
xlabel('Number of iterations')
ylabel('Dual step norm')
xlim([3 length(Fid)])


alphaArr = reshape(DP*u(1:m*n),[m,n]);
betaArr = reshape(DP*u(1+m*n:end),[m,n]);

estACtv = alphaArr*NptodB/100;
estBAtv = 2*(betaArr-1);

fprintf('AC: %.2f +/- %.2f\n', mean(estACtv(:),'omitnan'), ...
    std(estACtv(:), [] ,'omitnan'));
fprintf('B/A: %.2f +/- %.2f\n', mean(estBAtv(:),'omitnan'), ...
    std(estBAtv(:), [] ,'omitnan'));


figure; imagesc(xP*1e3,zP*1e3,estACtv); colorbar;
clim(attRange);
title('\alpha_0 in \alpha(f) = \alpha_0 \times f^2 dB/cm')
axis image
colormap turbo; colorbar;
xlabel('Lateral distance (mm)');
ylabel('Depth (mm)');
% hold on
% rectangle('Position',[0-radiusDisk,centerDepth-radiusDisk,...
% 2*radiusDisk,2*radiusDisk]*1000, 'Curvature',1,...
% 'EdgeColor','b', 'LineStyle','--', 'LineWidth',2)
% hold off


figure; imagesc(xP*1e3,zP*1e3,estBAtv); colorbar;
clim(baRange);
title('B/A');
axis image
colormap pink; colorbar;
xlabel('Lateral distance (mm)');
ylabel('Depth (mm)');
% hold on
% rectangle('Position',[0-radiusDisk,centerDepth-radiusDisk,...
% 2*radiusDisk,2*radiusDisk]*1000, 'Curvature',1,...
% 'EdgeColor','b', 'LineStyle','--', 'LineWidth',2)
% hold off
pause(0.1)


% save_all_figures_to_directory(fullfile(baseDir,'24-09-19'),'fig')
% close all

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
