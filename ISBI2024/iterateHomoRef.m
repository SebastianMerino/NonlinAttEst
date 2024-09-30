% New methods with frequency compounding. Requires three transmissions
% Used for showing results

setup;
baseDir = 'C:\Users\sebas\Documents\Data\Nonlinearity\homoAttMismatch\fn2';
resultsDir = fullfile(baseDir,'results');
[~,~,~] = mkdir(resultsDir);

% Auxiliar variables
NptodB = 20*log10(exp(1));
radiusDisk = (9)*1e-3;
centerDepth = 22.5e-3;

baRange = [5 13];
attRange = [0.09,0.17];


alphaRefVec = 0.08:0.02:0.22;
iSim = 1;
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
blockParams.blockSize = [20 20]*wl;
blockParams.overlap = 0.8;
blockParams.zlim = [0.8; 5.5]/100;
blockParams.xlim = [-2.5; 2.5]/100;

freqVec = [4,5,6]; % FRECUENCIES FOR FILTERING

%% For loop
for iSim=1:length(alphaRefVec)
    %%
    alphaRef = alphaRefVec(iSim);
    alphaInit = alphaRef;
    medium.alphaR = alphaRef/NptodB*100; % alpha0 in dB/100/MHz2

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
    end

    % Initialization
    [Xmesh,Zmesh] = meshgrid(xP,zP);
    alphaL = ones(size(Xmesh))*alphaInit;
    betaL = ones(size(Xmesh))*(1+baInit/2);
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


    %% Local maps with regularization
    muLocal = 0.1;
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

    metricsGN(iSim) = getMetrics(estACinst,estBAinst,'GN',alphaRef);

    %% ADMM
    % Optimizes F(u) + R(v)
    % Hyperparameters
    [m,n,p] = size(bzf);
    tol = 1e-3;
    muAlpha = 0.01; muBeta = 0.01;
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

    figure; imagesc(xP*1e3,zP*1e3,estBAtv); colorbar;
    clim(baRange);
    title("B/A = "+num2str(baMean,2)+"+/-"+num2str(baStd,2));
    axis image
    colormap pink; colorbar;
    xlabel('Lateral distance (mm)');
    ylabel('Depth (mm)');

    figure;
    im = imagesc(xP*1e3,zP*1e3,estACtv); colorbar;
    clim(attRange);
    title("\alpha_0 = "+num2str(acMean,2)+"+/-"+num2str(acStd,2));
    axis image
    colormap turbo; colorbar;
    xlabel('Lateral distance (mm)');
    ylabel('Depth (mm)');
    pause(0.1)

    metricsADMM(iSim) = getMetrics(estACtv,estBAtv,'ADMM',alphaRef);

    %%
    save_all_figures_to_directory(resultsDir,...
        char("sim"+iSim+"fig"))
    close all

end
%%
T = [struct2table(metricsGN);struct2table(metricsADMM)];
writetable(T,fullfile(resultsDir,'table.xlsx'))

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