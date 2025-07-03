% New methods with frequency compounding. Requires three transmissions
% Used for showing iterating and showing results for ba6inc12

setup;
baseDir = "Q:\smerino\Nonlinearity\AC_UiX_new\bf";
resultsDir = "Q:\smerino\Nonlinearity\resultsJASA\newSimulation\ba6inc12Ref2p0";
[~,~,~] = mkdir(resultsDir);
refDir = "Q:\smerino\Nonlinearity\AC_UiX_new\bf";

% Auxiliar variables
NptodB = 20*log10(exp(1));
radiusDisk = (9)*1e-3;
centerDepth = 22.5e-3;

imPosition = [100 200 200 250];
baRange = [4 13];
attRange = [0.08,0.22];

alphaIncVec = 8:2:20;

alphaInit = 0.1;
baInit = 6;

%% Preparing data
% Known variables
medium.v = 5;
medium.betaR = 1 + 6/2;
medium.alphaR = 0.1/NptodB*100; % alpha0 in dB/100/MHz2

% Filtering parameters
filterParams.freqTol = 0.5;
filterParams.nCycles = 10; % Number of cycles of the initial filter

% Subsampling parameters
wl = 1540/6e6; % Mean central frequency
blockParams.blockSize = [25,25]*wl;
blockParams.overlap = 0.8;
blockParams.zlim = [0.5; 5.4]/100;
blockParams.xlim = [-2.5; 2.5]/100;
blockParams.downFactor = 20;
freqVec = [5,6,7]; % FRECUENCIES FOR FILTERING

iSim = 7;
alphaInc = alphaIncVec(iSim);

%% For loop
for iSim=1:length(alphaIncVec)
    alphaInc = alphaIncVec(iSim);
    % ---------------------------------------------------------------------- %
    % ---------------------------------------------------------------------- %
    % ---------------------------------------------------------------------- %
    %% Measurements, IUS version
    freq = 5;
    alphaStr = num2str(alphaInc,"%02d");
    fileSam = "RFfn2_PWNE"+freq+"MHz_sam_att0p1inc0p"+alphaStr+ ...
            "f20_BA6inc12_nc10_400kPa";
    fileRef = "RFfn2_PWNE"+freq+"MHz_ref_att0p1f20_BA6_nc10_400kPa";
    
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
    filterParams.freqC = freq;
    [bz,zbz,xB,zB] = getMeasurementsIUS(medium,filterParams,blockParams);
    X = permute(zbz,[3 1 2]);
    Y = permute(bz,[3 1 2]);
    X = X(:); Y = Y(:);

    %% Gauss-Newton with LM
    [m,n,p] = size(bz);
    
    tol = 1e-3;
    maxIte = 200;
    muAlpha = 3; muBeta = 3;
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
    exTime = toc;
    alphaArr = theta(1:n*m);
    betaArr = theta(n*m+1:end);
    
    estAClm = reshape(alphaArr/freq^2 /100*NptodB,[m,n]);
    estBAlm = reshape(2*(betaArr-1),[m,n]);

    %% Local maps with regularization
    muLocal = 0.001;
    dzP = zB(2)-zB(1);
    izP = round(zB./dzP);
    factorq = izP(1)./izP;
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

    %% ROI and metrics
    %  Masks
    xBm = medium.x; zBm = medium.z;
    [Xq, Zq] = meshgrid(xBm,zBm);
    Lx = radiusDisk*1.1; Lz = radiusDisk*1.4; 
    cx = radiusDisk*1.5;

    inc = maskRect(xBm, zBm, 0, centerDepth, Lx, Lz);
    back = maskRect(xBm, zBm, -cx, centerDepth, Lx/2, Lz) | ...
        maskRect(xBm, zBm, cx, centerDepth, Lx/2, Lz);

    % Interp
    [Xmesh,Zmesh] = meshgrid(xB,zB(2:end));
    baInterp = interp2(Xmesh,Zmesh,estBAinst,Xq,Zq);
    acInterp = interp2(Xmesh,Zmesh,estACinst,Xq,Zq);

    % Metrics
    metricsGN(iSim) = getMetrics(acInterp,baInterp,inc,back,'IUS', ...
        alphaInc/100,exTime,ite);

    %% Ideal map
    incBm = Xq.^2 + (Zq-centerDepth).^2 < radiusDisk^2;
    idealAC = 0.1*ones(size(Xq));
    idealAC(incBm) = alphaInc/100;
    idealBA = 6*ones(size(Xq));
    idealBA(incBm) = 12;

    figure('Position',imPosition); 
    im = imagesc(xBm*1e2,zBm*1e2,idealAC); colorbar;
    clim(attRange);
    title('\alpha_0 in \alpha(f) = \alpha_0 \times f^2')
    axis equal
    xlim([xB(1),xB(end)]*100)
    ylim([zB(1),zB(end)]*100)
    colormap turbo; colorbar;
    xlabel('Lateral [cm]');
    ylabel('Depth [cm]');

    figure('Position',imPosition);  
    imagesc(xBm*1e2,zBm*1e2,idealBA); colorbar;
    clim(baRange);
    title('B/A');
    axis equal
    xlim([xB(1),xB(end)]*100)
    ylim([zB(1),zB(end)]*100)
    colormap pink; colorbar;
    xlabel('Lateral [cm]');
    ylabel('Depth [cm]');
    pause(0.1)


    %% Plots
    figure('Position',imPosition); 
    im = imagesc(xB*1e2,zB*1e2,estACinst); colorbar;
    clim(attRange);
    title('\alpha_0 in \alpha(f) = \alpha_0 \times f^2')
    axis image
    colormap turbo; colorbar;
    xlabel('Lateral [cm]');
    ylabel('Depth [cm]');
    hold on
    rectangle('Position',[0-radiusDisk,centerDepth-radiusDisk,...
        2*radiusDisk,2*radiusDisk]*1e2, 'Curvature',1,...
        'EdgeColor','w', 'LineStyle','--', 'LineWidth',2)
    rectangle('Position',[-Lx/2,centerDepth-Lz/2,...
        Lx,Lz]*1e2, 'EdgeColor','k', 'LineStyle','--', 'LineWidth',2)
    rectangle('Position',[cx-Lx/4,centerDepth-Lz/2,...
        Lx/2,Lz]*1e2, 'EdgeColor','k', 'LineStyle','--', 'LineWidth',2)
    rectangle('Position',[-cx-Lx/4,centerDepth-Lz/2,...
        Lx/2,Lz]*1e2, 'EdgeColor','k', 'LineStyle','--', 'LineWidth',2)
    hold off

    figure('Position',imPosition);  
    imagesc(xB*1e2,zB(2:end)*1e2,estBAinst); colorbar;
    clim(baRange);
    title('B/A');
    axis image
    colormap pink; colorbar;
    xlabel('Lateral [cm]');
    ylabel('Depth [cm]');
    hold on
    rectangle('Position',[0-radiusDisk,centerDepth-radiusDisk,...
        2*radiusDisk,2*radiusDisk]*1e2, 'Curvature',1,...
        'EdgeColor','w', 'LineStyle','--', 'LineWidth',2)
    rectangle('Position',[-Lx/2,centerDepth-Lz/2,...
        Lx,Lz]*1e2, 'EdgeColor','k', 'LineStyle','--', 'LineWidth',2)
    rectangle('Position',[cx-Lx/4,centerDepth-Lz/2,...
        Lx/2,Lz]*1e2, 'EdgeColor','k', 'LineStyle','--', 'LineWidth',2)
    rectangle('Position',[-cx-Lx/4,centerDepth-Lz/2,...
        Lx/2,Lz]*1e2, 'EdgeColor','k', 'LineStyle','--', 'LineWidth',2)
    hold off
    pause(0.1)

    % ---------------------------------------------------------------------- %
    % ---------------------------------------------------------------------- %
    % ---------------------------------------------------------------------- %
    %% Getting B/A
    bzf = [];
    for iFreq = 1:length(freqVec)
        freq = freqVec(iFreq);
        alphaStr = num2str(alphaInc,"%02d");
        fileSam = "RFfn2_PWNE"+freq+"MHz_sam_att0p1inc0p"+alphaStr+ ...
                "f20_BA6inc12_nc10_400kPa";
        fileRef = "RFfn2_PWNE"+freq+"MHz_ref_att0p1f20_BA6_nc10_400kPa";

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

        % idz = medium.z>0.5e-2&medium.z<5.5e-2; % true(size(medium.z)); %
        % bmode = db(hilbert(medium.rfH(idz,:,1)));
        % bmode = bmode - max(bmode (:));
        % figure,
        % imagesc(medium.x*100,medium.z(idz)*100,bmode, [-40 0])
        % colormap gray
        % axis image
        % colorbar
        % title('B-mode')
    end

    %% Initialization
    [m,n,p] = size(bzf);
    alphaL = ones(m,n)*alphaInit;
    betaL = ones(m,n)*(1+baInit/2);
    u0 = [alphaL(:)*100/NptodB;betaL(:)];

    %% ADMM
    % Optimizes F(u) + R(v)
    % Hyperparameters
    [m,n,p] = size(bzf);
    tol = 1e-4;
    muAlpha = 10^(-2); muBeta = 10^(-2.5);
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
    exTime = toc;

    alphaArr = reshape(DP*u(1:m*n),[m,n]);
    betaArr = reshape(DP*u(1+m*n:end),[m,n]);

    estACtv = alphaArr*NptodB/100;
    estBAtv = 2*(betaArr-1);
    %% ROI and metrics
    %  Masks
    xBm = medium.x; zBm = medium.z;
    [Xq, Zq] = meshgrid(xBm,zBm);
    inc = maskRect(xBm, zBm, 0, centerDepth, Lx, Lz);
    back = maskRect(xBm, zBm, -cx, centerDepth, Lx/2, Lz) | ...
        maskRect(xBm, zBm, cx, centerDepth, Lx/2, Lz);

    % Interp
    [Xmesh,Zmesh] = meshgrid(xP,zP);
    baInterp = interp2(Xmesh,Zmesh,estBAtv,Xq,Zq);
    acInterp = interp2(Xmesh,Zmesh,estACtv,Xq,Zq);

    % Metrics
    metricsADMM(iSim) = getMetrics(acInterp,baInterp,inc,back,'ADMM', ...
        alphaInc/100,exTime,ite);

    %% Plots
    figure('Position',imPosition); 
    im = imagesc(xP*1e2,zP*1e2,estACtv); colorbar;
    clim(attRange);
    title('\alpha_0 in \alpha(f) = \alpha_0 \times f^2')
    axis image
    colormap turbo; colorbar;
    xlabel('Lateral [cm]');
    ylabel('Depth [cm]');
    hold on
    rectangle('Position',[0-radiusDisk,centerDepth-radiusDisk,...
        2*radiusDisk,2*radiusDisk]*1e2, 'Curvature',1,...
        'EdgeColor','w', 'LineStyle','--', 'LineWidth',2)
    rectangle('Position',[-Lx/2,centerDepth-Lz/2,...
        Lx,Lz]*1e2, 'EdgeColor','k', 'LineStyle','--', 'LineWidth',2)
    rectangle('Position',[cx-Lx/4,centerDepth-Lz/2,...
        Lx/2,Lz]*1e2, 'EdgeColor','k', 'LineStyle','--', 'LineWidth',2)
    rectangle('Position',[-cx-Lx/4,centerDepth-Lz/2,...
        Lx/2,Lz]*1e2, 'EdgeColor','k', 'LineStyle','--', 'LineWidth',2)
    hold off

    figure('Position',imPosition);  
    imagesc(xP*1e2,zP(2:end)*1e2,estBAtv); colorbar;
    clim(baRange);
    title('B/A');
    axis image
    colormap pink; colorbar;
    xlabel('Lateral [cm]');
    ylabel('Depth [cm]');
    hold on
    rectangle('Position',[0-radiusDisk,centerDepth-radiusDisk,...
        2*radiusDisk,2*radiusDisk]*1e2, 'Curvature',1,...
        'EdgeColor','w', 'LineStyle','--', 'LineWidth',2)
    rectangle('Position',[-Lx/2,centerDepth-Lz/2,...
        Lx,Lz]*1e2, 'EdgeColor','k', 'LineStyle','--', 'LineWidth',2)
    rectangle('Position',[cx-Lx/4,centerDepth-Lz/2,...
        Lx/2,Lz]*1e2, 'EdgeColor','k', 'LineStyle','--', 'LineWidth',2)
    rectangle('Position',[-cx-Lx/4,centerDepth-Lz/2,...
        Lx/2,Lz]*1e2, 'EdgeColor','k', 'LineStyle','--', 'LineWidth',2)
    hold off
    pause(0.1)

    %%
    save_all_figures_to_directory(resultsDir,...
        char("sim"+iSim+"fig"),'svg')
    close all
end

%%
T = [struct2table(metricsGN);struct2table(metricsADMM)];
writetable(T,fullfile(resultsDir,'table.xlsx'))

%% Utility functions
function metrics = getMetrics(AC,BA,inc,back,method,alphaInc,time,ite)
metrics.AcIncMean = mean(AC(inc), 'omitnan');
metrics.AcIncStd = std(AC(inc), [], 'omitnan');
metrics.AcBackMean = mean(AC(back), 'omitnan');
metrics.AcBackStd = std(AC(back), [], 'omitnan');
metrics.BaIncMean = mean(BA(inc), 'omitnan');
metrics.BaIncStd = std(BA(inc), [], 'omitnan');
metrics.BaBackMean = mean(BA(back), 'omitnan');
metrics.BaBackStd = std(BA(back), [], 'omitnan');
metrics.method = method;
metrics.alphaInc = alphaInc;
metrics.execTime = time;
metrics.nIter = ite;
end