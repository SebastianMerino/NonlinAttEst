% New methods with frequency compounding. Requires three transmissions
% Used for showing results

startup;
baseDir = 'Q:\smerino\Nonlinearity\AC_UiX_new\homog\bf';
resultsDir = 'Q:\smerino\Nonlinearity\resultsJASA\newSimulation\homoBA9';
[~,~,~] = mkdir(resultsDir);

% Auxiliar variables
NptodB = 20*log10(exp(1));
radiusDisk = (9)*1e-3;
centerDepth = 22.5e-3;

baRange = [5 13];
attRange = [0.05,0.15];
imPosition = [100 200 200 250];

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
wl = 1540/6e6; % Mean central frequency
blockParams.blockSize = [25 25]*wl;
blockParams.overlap = 0.8;
blockParams.zlim = [1; 5]/100;
blockParams.xlim = [-2.5; 2.5]/100;
blockParams.downFactor = 20;

freqVec = [5,6,7]; % FRECUENCIES FOR FILTERING
iSim = 7;

%% For loop
for iSim=1:length(alphaRefVec)
    %%
    zMax = -(iSim-1)/7*3+5.5;
    blockParams.zlim(2) = zMax/100;
    alphaRef = alphaRefVec(iSim);
    alphaInit = alphaRef;
    medium.alphaR = alphaRef/NptodB*100; % alpha0 in dB/100/MHz2

    %% Measurements, IUS version
    freq = 5;
    fileSam = "RFfn2_PWNE"+freq+"MHz_sam_att0p10f20_BA9_nc10_400kPa";
    fileRef = "RFfn2_PWNE"+freq+"MHz_ref_att0p"+...
        num2str(round(100*alphaRef), '%02d')+"f20_BA6_nc10_400kPa";
    
    % Sample
    sample = load(fullfile(baseDir,fileSam));
    medium.z = sample.z';
    medium.x = sample.x;
    medium.fs = sample.fs;
    medium.rfL = sample.rf1(:,:,1:2);
    medium.rfH = sample.rf2(:,:,1:2);
    clear sample
    
    % Reference
    ref = load(fullfile(baseDir,fileRef));
    medium.rfLR = ref.rf1(:,:,1:2);
    medium.rfHR = ref.rf2(:,:,1:2);
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

    baMean = mean(estBAlm(:),'all','omitnan');
    baStd = std(estBAlm(:),[],'all','omitnan');
    acMean = mean(estAClm(:),'all','omitnan');
    acStd = std(estAClm(:),[],'all','omitnan');

    figure('Position',imPosition); 
    imagesc(xB*1e2,zB*1e2,estBAlm); colorbar;
    clim(baRange);
    title("B/A = "+sprintf("%.1f",baMean)+"\pm"+sprintf("%.1f",baStd));
    axis image
    colormap pink; colorbar;
    xlabel('Lateral [cm]');
    ylabel('Depth [cm]');

    figure('Position',imPosition); 
    imagesc(xB*1e2,zB*1e2,estAClm); colorbar;
    clim(attRange);
    title("\alpha_0 = "+sprintf("%.2f",acMean)+"\pm"+sprintf("%.2f",acStd));
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
        fileSam = "RFfn2_PWNE"+freq+"MHz_sam_att0p10f20_BA9_nc10_400kPa";
        fileRef = "RFfn2_PWNE"+freq+"MHz_ref_att0p"+...
            num2str(round(100*alphaRef), '%02d')+"f20_BA6_nc10_400kPa";

        % Sample
        sample = load(fullfile(baseDir,fileSam));
        medium.z = sample.z';
        medium.x = sample.x;
        medium.fs = sample.fs;
        medium.rfL = sample.rf1(:,:,1:2);
        medium.rfH = sample.rf2(:,:,1:2);
        clear sample

        % Reference
        ref = load(fullfile(baseDir,fileRef));
        medium.rfLR = ref.rf1(:,:,1:2);
        medium.rfHR = ref.rf2(:,:,1:2);
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
    muAlpha = 10^(-1); muBeta = 10^(-1.5);
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

    baMean = mean(estBAtv(:),'all','omitnan');
    baStd = std(estBAtv(:),[],'all','omitnan');
    acMean = mean(estACtv(:),'all','omitnan');
    acStd = std(estACtv(:),[],'all','omitnan');

    figure('Position',imPosition); 
    imagesc(xP*1e2,zP*1e2,estBAtv); colorbar;
    clim(baRange);
    title("B/A = "+sprintf("%.1f",baMean)+"\pm"+sprintf("%.1f",baStd));
    axis image
    colormap pink; colorbar;
    xlabel('Lateral [cm]');
    ylabel('Depth [cm]');

    figure('Position',imPosition); 
    im = imagesc(xP*1e2,zP*1e2,estACtv); colorbar;
    clim(attRange);
    title("\alpha_0 = "+sprintf("%.2f",acMean)+"\pm"+sprintf("%.2f",acStd));
    axis image
    colormap turbo; colorbar;
    xlabel('Lateral [cm]');
    ylabel('Depth [cm]');
    pause(0.1)

    metricsADMM(iSim) = getMetrics(estACtv,estBAtv,'ADMM',alphaRef);

    %%
    save_all_figures_to_directory(resultsDir,...
        char("sim"+iSim+"fig"),'svg')
    close all

end
%%
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