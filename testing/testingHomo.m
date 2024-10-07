% New methods with frequency compounding. Requires three transmissions
% Used for showing results

setup;
baseDir = 'C:\Users\sebas\Documents\Data\Nonlinearity\homoAttMismatch\fn2';
% Auxiliar variables
NptodB = 20*log10(exp(1));
radiusDisk = (9)*1e-3;
centerDepth = 22.5e-3;

baRange = [5 13];
attRange = [0.09,0.17];

alphaRef = 0.14;

alphaInit = alphaRef;
betaInit = 6;

%% Preparing data
% Known variables
medium.v = 5;
medium.betaR = 1 + 6/2;             
medium.alphaR = alphaRef/NptodB*100; % alpha0 in dB/100/MHz2

% Filtering parameters
filterParams.freqC = 5;
filterParams.freqTol = 0.5;
filterParams.nCycles = 10; % Number of cycles of the initial filter

% Subsampling parameters
wl = 1540/5e6; % Mean central frequency
blockParams.blockSize = [25 25]*wl; 
blockParams.overlap = 0.8;
blockParams.zlim = [0.8; 5.5]/100;
blockParams.xlim = [-2.5; 2.5]/100;

freqVec = [4,5,6]; % FRECUENCIES FOR FILTERING

%% Getting B/A
bzf = [];
for iFreq = 1:length(freqVec)
    freq = freqVec(iFreq);
    fileSam = "RFfn2_PWNE"+freq+"MHz_samBA9_att0p10f2_nc10_400kPa";
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

% Initialization
[Xmesh,Zmesh] = meshgrid(xP,zP);
alphaL = ones(size(Xmesh))*alphaInit;
betaL = ones(size(Xmesh))*(1+betaInit/2);
u0 = [alphaL(:)*100/NptodB;betaL(:)];


%% Gauss-Newton with LM
[m,n,p] = size(bzf);
tol = 1e-3;
maxIte = 200;
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

%% Local maps with regularization
muLocal = 0.1; %0.001;
dzP = zP(2)-zP(1);
izP = round(zP./dzP);
factorq = izP(1)./izP ;
P = sparse(tril(ones(m-1)));
P = P./izP(2:end);
P = kron(speye(n),P);

estBAcum = estBAlm - estBAlm(1,:).*factorq; 
estBAcum = estBAcum(2:end,:);
estBAinst = IRLS_TV(estBAcum(:),P,muLocal,m-1,n,tol,[],ones((m-1)*n,1));
estBAinst = reshape(estBAinst,m-1,n);

estACcum = estAClm - estAClm(1,:).*factorq; 
estACcum = estACcum(2:end,:);
estACinst = IRLS_TV(estACcum(:),P,muLocal,m-1,n,tol,[],ones((m-1)*n,1));
estACinst = reshape(estACinst,m-1,n);

baMean = mean(estBAinst(:),'omitnan');
baStd = std(estBAinst(:),[],'omitnan');
acMean = mean(estACinst(:),'omitnan');
acStd = std(estACinst(:),[],'omitnan');

figure; imagesc(xP*1e3,zP(2:end)*1e3,estBAinst); colorbar;
clim(baRange);
title("B/A = "+num2str(baMean,2)+"+/-"+num2str(baStd,2));
axis image
colormap pink; colorbar;
xlabel('Lateral distance (mm)');
ylabel('Depth (mm)');

figure; 
im = imagesc(xP*1e3,zP*1e3,estACinst); colorbar;
clim(attRange);
title("\alpha_0 = "+num2str(acMean,2)+"+/-"+num2str(acStd,2));
axis image
colormap turbo; colorbar;
xlabel('Lateral distance (mm)');
ylabel('Depth (mm)');
pause(0.1)

%% ADMM
% Optimizes F(u) + R(v)
% Hyperparameters
[m,n,p] = size(bzf);
tol = 1e-3;
muAlpha = 0.01; muBeta = 0.01;
% muAlpha = 0.1; muBeta = 0.01;
rho = 0.1;
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
figure; tiledlayout(1,2)
t1 = nexttile;         t2 = nexttile;
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

    if mod(ite,10) ==1 && ite<50
        alphaArr = reshape(DP*u(1:m*n),[m,n]);
        betaArr = reshape(DP*u(1+m*n:end),[m,n]);
        estACtv = alphaArr*NptodB/100;
        estBAtv = 2*(betaArr-1);

        set(gcf, 'currentaxes', t1);
        imagesc(t1,xP*1e3,zP*1e3,estACtv); colorbar;
        clim(attRange);
        title('\alpha_0 in \alpha(f) = \alpha_0 \times f^2 dB/cm')
        axis image
        colorbar(t1);
        xlabel('Lateral distance (mm)');
        ylabel('Depth (mm)');

        set(gcf, 'currentaxes', t2);
        imagesc(t2,xP*1e3,zP*1e3,estBAtv); colorbar;
        clim(baRange);
        title('B/A');
        axis image
        colormap(t2,pink); colorbar;
        xlabel('Lateral distance (mm)');
        ylabel('Depth (mm)');

        colormap(t1,turbo); 
        pause(0.1)

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


[Xmesh,Zmesh] = meshgrid(xP,zP);
inc = Xmesh.^2 + (Zmesh-centerDepth).^2 < (radiusDisk-1e-3).^2;
back = Xmesh.^2 + (Zmesh-centerDepth).^2 > (radiusDisk+1e-3).^2;

fprintf('AC inc: %.2f +/- %.2f\n', mean(estACtv(inc),'omitnan'), ...
std(estACtv(inc), [] ,'omitnan'));
fprintf('AC back: %.2f +/- %.2f\n', mean(estACtv(back),'omitnan'), ...
std(estACtv(back), [] ,'omitnan'));

fprintf('B/A inc: %.2f +/- %.2f\n', mean(estBAtv(inc),'omitnan'), ...
std(estBAtv(:), [] ,'omitnan'));
fprintf('B/A back: %.2f +/- %.2f\n', mean(estBAtv(back),'omitnan'), ...
std(estBAtv(:), [] ,'omitnan'));


figure; imagesc(xP*1e3,zP*1e3,estACtv); colorbar;
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


figure; imagesc(xP*1e3,zP*1e3,estBAtv); colorbar;
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


% save_all_figures_to_directory(fullfile(baseDir,'24-09-19'),'fig')
% close all

% ---------------------------------------------------------------------- %
% ---------------------------------------------------------------------- %