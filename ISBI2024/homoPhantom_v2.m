% New methods with frequency compounding. Requires three transmissions
% Used for showing results. Attenuation gamma can be set

setup;
baseDir = 'C:\Users\sebas\Documents\Data\Nonlinearity\01-Oct-24\bf';

% Auxiliar variables
NptodB = 20*log10(exp(1));
radiusDisk = (9)*1e-3;
centerDepth = 22.5e-3;

baRange = [4 12];
attRange = [0.05,0.2];

alphaInit = 0.134; % change
baInit = 10.5;

c0 = 1470;
%% Preparing data
freqVec = [5,6,7]; % FRECUENCIES FOR FILTERING
gammaAtt = 1.6;

% Known variables
medium.v = 4; % 3.3829 in fit, 4 visually

% % Reference adriana
medium.betaR = 1 + 5.4/2;
medium.alphaR = 0.3504/NptodB*100; %alpha0 in dB/m/MHz^2
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
blockParams.blockSize = [25 25]*wl;
blockParams.overlap = 0.8;
blockParams.zlim = [2;5.5]/100;
% blockParams.zlim = [0.5;1]/100;
blockParams.xlim = [-1.5; 1.5]/100;

vArr = [5,4.5,4.4];
%% Getting B/A
bzf = [];
for iFreq = 1:length(freqVec)
    %%
    medium.v = vArr(iFreq);
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
    [bz,xP,zP] = getMeasurements_v2(medium,filterParams,blockParams,medium.alphaRpower);
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

    ejeX = medium.rfLR(idz,:,:); ejeY = medium.rfHR(idz,:,:);
    [values,c] = hist3([ejeX(:) ejeY(:)],[101 101]);
    figure,
    imagesc(c{1},c{2},db(values'))
    axis xy
    hold on, plot(-2e4:2e4,(-2e4:2e4)*medium.v,'k-'), hold off
    colorbar

end

%% Measurements
% bLim = [2,10];
bLim = [2,5.5];

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
%%
iz = 10;
betaIdeal = 1+10.5/2;
alphaIdeal = 0.3287;
% betaIdeal = 1+4.8/2;
% alphaIdeal = 0.11;
freqInterp = linspace(4,8,100);
alphaF = alphaIdeal/NptodB*100*freqInterp.^gammaAtt;
bzIdeal = betaIdeal*( 1-exp(-2*alphaF*zP(iz)) )./(alphaF*zP(iz));

bzf5 = bzf(iz,:,1);
bzf6 = bzf(iz,:,2);
bzf7 = bzf(iz,:,3);
figure,
plot(freqVec,[bzf5;bzf6;bzf7])
hold on
plot(freqInterp,bzIdeal, 'k--', 'LineWidth',2)
hold off
title(sprintf('Measurements at z = %.2f cm',zP(iz)*100))

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
loss(1) = 0.5*norm(modelFreq(u,zP,freqVec,gammaAtt) - bzf(:))^2;
ite = 1;
tic
while true
    jcb = jacobianFreq(u,zP,freqVec,gammaAtt);
    res = modelFreq(u,zP,freqVec,gammaAtt) - bzf(:);
    [step,~] = pcg(jcb'*jcb + regMatrix,jcb'*-res, tol, 200);
    % step = (jcb'*jcb + regMatrix)\(jcb'*-res);
    u = u + step;

    loss(ite+1) = 0.5*norm(modelFreq(u,zP,freqVec,gammaAtt) - bzf(:))^2;
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
title('\alpha_0 in \alpha(f) = \alpha_0 \times f^\gamma')
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
Fid(1) = 1/2*norm( modelFreq(u,zP,freqVec, gammaAtt) - bzf(:) )^2;
Reg(1) = muAlpha*TVcalc_isotropic(DP*u(1:m*n),m,n,mask) + ...
    muBeta*TVcalc_isotropic(DP*u(m*n+1:end),m,n,mask);
Dual(1) = 0;
ite  = 0;
error = 1;
f = figure; tiledlayout(1,2)
t1 = nexttile;         t2 = nexttile;
tic
while ite < maxIte
    ite = ite + 1;

    % Fidelity step
    u = optimNonLinearGNFreq(zP,freqVec,bzf(:), speye(2*m*n),v-w, rho, ...
        tol,u, gammaAtt);

    % Regularization step
    v = optimAdmmTvDy(Ii,Id,u+w, muAlpha/rho,muBeta/rho ,m,n,tol,mask,DP);

    % Dual update
    w = w + u - v;

    % Loss
    Fid(ite+1) = 1/2*norm( modelFreq(u,zP,freqVec, gammaAtt) - bzf(:) )^2;
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

        set(f, 'currentaxes', t1);
        imagesc(t1,xP*1e3,zP*1e3,estACtv); colorbar;
        clim(attRange);
        title('\alpha_0 in \alpha(f) = \alpha_0 \times f^2 dB/cm')
        axis image
        colorbar(t1);
        xlabel('Lateral distance (mm)');
        ylabel('Depth (mm)');

        set(f, 'currentaxes', t2);
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


% [Xmesh,Zmesh] = meshgrid(xP,zP);
% inc = Xmesh.^2 + (Zmesh-centerDepth).^2 < (radiusDisk-1e-3).^2;
% back = Xmesh.^2 + (Zmesh-centerDepth).^2 > (radiusDisk+1e-3).^2;
%
% fprintf('AC inc: %.2f +/- %.2f\n', mean(estACtv(inc),'omitnan'), ...
% std(estACtv(inc), [] ,'omitnan'));
% fprintf('AC back: %.2f +/- %.2f\n', mean(estACtv(back),'omitnan'), ...
% std(estACtv(back), [] ,'omitnan'));
% fprintf('B/A inc: %.2f +/- %.2f\n', mean(estBAtv(inc),'omitnan'), ...
% std(estBAtv(:), [] ,'omitnan'));
% fprintf('B/A back: %.2f +/- %.2f\n', mean(estBAtv(back),'omitnan'), ...
% std(estBAtv(:), [] ,'omitnan'));
fprintf('AC: %.3f +/- %.3f\n', mean(estACtv(:),'omitnan'), ...
    std(estACtv(:), [] ,'omitnan'));
fprintf('B/A: %.2f +/- %.2f\n', mean(estBAtv(:),'omitnan'), ...
    std(estBAtv(:), [] ,'omitnan'));


figure; imagesc(xP*1e3,zP*1e3,estACtv); colorbar;
clim(attRange);
title('\alpha_0 in \alpha(f) = \alpha_0 \times f^\gamma dB/cm')
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