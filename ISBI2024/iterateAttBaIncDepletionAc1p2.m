% New methods with frequency compounding. Requires three transmissions
% Used for showing results

setup;
baseDir = "Q:\smerino\Nonlinearity\attIncNonQuadratic";
resultsDir = "Q:\smerino\Nonlinearity\resultsJASA\ba6inc12ac1p2Depletionref2p0";
[~,~,~] = mkdir(resultsDir);
refDir = "Q:\smerino\Nonlinearity\newRef";

% Auxiliar variables
NptodB = 20*log10(exp(1));
radiusDisk = (9)*1e-3;
centerDepth = 22.5e-3;

imPosition = [100 200 250 300];
baRange = [4 13];
attRange = [0.08,0.22];

alphaIncVec = 8:2:14;

alphaInit = 0.1;
baInit = 6;

%% Preparing data
% Known variables
medium.v = 5;
medium.betaR = 1 + 6/2;
medium.alphaRcoeff = 0.1/NptodB*100; % alpha0 in dB/m
medium.alphaRpower = 2.0;

% Filtering parameters
filterParams.freqC = 5;
filterParams.freqTol = 0.5;
filterParams.nCycles = 10; % Number of cycles of the initial filter

% Subsampling parameters
wl = 1540/5e6; % Mean central frequency
blockParams.blockSize = [15 15]*wl;
blockParams.overlap = 0.8;
blockParams.zlim = [0.5; 5.5]/100;
blockParams.xlim = [-2.5; 2.5]/100;
blockParams.downFactor = 20;
freqVec = [4,5,6]; % FRECUENCIES FOR FILTERING

% Attenuation estimation parameters
blockParamsAcs.xInf = -2.5;
blockParamsAcs.xSup = 2.5;
blockParamsAcs.zInf = 0.5;
blockParamsAcs.zSup = 5.5;
blockParamsAcs.blocksize = [15 15]*wl;
blockParamsAcs.overlap = 0.8;
freqL = 4e6; freqH = 6e6;
tol = 1e-3;
muRsld = 1e2;

iSim = 4;
alphaInc = alphaIncVec(iSim);

%% For loop
for iSim = 1:length(alphaIncVec)
    alphaInc = alphaIncVec(iSim);
    %% ================================================================ %%
    %% Measurements, DEPLETION
    bzf = [];

    iFreq = 2;
    freq = freqVec(iFreq);
    alphaStr = num2str(alphaInc,"%02d");
    fileSam = "RFfn2_PWNE"+freq+"MHz_sam_att0p1inc0p"+alphaStr+ ...
        "f12_BA6inc12_nc10_400kPa";
    fileRef = "RFfn2_PWNE"+freq+"MHz_ref_att0p1f20_BA6_nc10_400kPa";

    %% Sample
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

    %% Building ideal attenuation maps from sample
    % Alpha coefficients in Np/m
    alphaCoeffBack = 0.1/NptodB*100; 
    alphaCoeffInc = alphaInc/NptodB;
    alphaPower = 1.2;

    attSamBack = alphaCoeffBack*(filterParams.freqC)^alphaPower;
    attSamInc = alphaCoeffInc*(filterParams.freqC)^alphaPower;

%    xPm = medium.x; zPm = medium.z;
    [X, Z] = meshgrid(xP,zP);
    attSam = attSamBack*ones(size(X));
    attSam(X.^2 + (Z-centerDepth).^2 < radiusDisk^2) = attSamInc;
    figure('Position',imPosition),
    imagesc(xP*100,zP*100,attSam)
    axis image
    colorbar
    title('Ideal Attenuation at 5 MHz')
    attSam = cumsum(attSam)./(1:length(zP))';
    estBAdep = (bz./( 1-exp(-2*attSam.*zP) ) .*attSam.*zP - 1)*2 ;


    %% Local maps with regularization
    [m,n] = size(bz);
    tol = 1e-3;
    muLocal = 0.01;
    dzP = zP(2)-zP(1);
    izP = round(zP./dzP);
    factorq = izP(1)./izP;
    P = sparse(tril(ones(m-1)));
    P = P./izP(2:end);
    P = kron(speye(n),P);

    estBAcum = estBAdep - estBAdep(1,:).*factorq;
    estBAcum = estBAcum(2:end,:);
    estBAinst = IRLS_TV(estBAcum(:),P,muLocal,m-1,n,tol,[],ones((m-1)*n,1));
    estBAinst = reshape(estBAinst,m-1,n);

    %% ROI and metrics
    %  Masks
    xPm = medium.x; zPm = medium.z;
    [Xq, Zq] = meshgrid(xPm,zPm);
    Lx = radiusDisk*1; Lz = radiusDisk*1.3;
    cx = radiusDisk*1.45;

    inc = maskRect(xPm, zPm, 0, centerDepth, Lx, Lz);
    back = maskRect(xPm, zPm, -cx, centerDepth, Lx/2, Lz) | ...
        maskRect(xPm, zPm, cx, centerDepth, Lx/2, Lz);

    % Interp
    [Xmesh,Zmesh] = meshgrid(xP,zP(2:end));
    baInterp = interp2(Xmesh,Zmesh,estBAinst,Xq,Zq);
    acInterp = zeros(size(baInterp));

    % Metrics
    metricsDep(iSim) = getMetrics(acInterp,baInterp,inc,back,'DepIdeal', ...
        alphaInc/100,0,1);

    %% Plots

    figure('Position',imPosition);
    imagesc(xP*1e2,zP(2:end)*1e2,estBAinst); colorbar;
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

    %% ================================================================ %%
    %% ATTENUATION ESTIMATION    
    % get Spectra
    [Sp,Sd,x_ACS,z_ACS,f] = getSpectrum(medium.rfL,medium.x*100,medium.z*100, ...
        medium.fs,blockParamsAcs);
    att_ref_map(1,1,:) = medium.alphaRcoeff/100*f.^medium.alphaRpower;
    [SpRef,SdRef,~,~,~] = getSpectrum(medium.rfLR,medium.x*100,medium.z*100, ...
        medium.fs,blockParamsAcs);

    % Setting up system
    L = (z_ACS(2) - z_ACS(1))/(1 - blockParamsAcs.overlap)/2;   % (cm)
    sld = log(Sp) - log(Sd);
    sldRef = log(SpRef) - log(SdRef);
    compensation = sldRef - 4*L*att_ref_map;
    range = f>freqL/1e6 & f<freqH/1e6;
    b = sld(:,:,range) - compensation(:,:,range);
    ufr = f(range);
    [m,n,p] = size(b);
    A1 = kron( 4*L*ufr , speye(m*n) );
    A2 = kron( ones(size(ufr)) , speye(m*n) );
    mask = ones(m,n,p);

    % Inversion
    tic
    [Bn,Cn] = AlterOpti_ADMM(A1,A2,b(:),muRsld,muRsld,m,n,tol,mask(:));
    toc
    BR = (reshape(Bn,m,n));
    CR = (reshape(Cn,m,n));

    % Plot and cumulative map
    attSam = (BR*filterParams.freqC + CR)*100; % Np/cm to Np/m
    figure('Position',imPosition),
    imagesc(xP*100,zP*100,attSam)
    axis image
    c = colorbar;
    c.Label.String = 'dB/cm';
    title('Attenuation from RSLD')
    colormap turbo
    attSam = cumsum(attSam)./(1:length(zP))';
    estBAdep = (bz./( 1-exp(-2*attSam.*zP) ) .*attSam.*zP - 1)*2 ;

    %% Local maps with regularization
    [m,n] = size(bz);
    tol = 1e-3;
    muLocal = 0.01;
    dzP = zP(2)-zP(1);
    izP = round(zP./dzP);
    factorq = izP(1)./izP;
    P = sparse(tril(ones(m-1)));
    P = P./izP(2:end);
    P = kron(speye(n),P);

    estBAcum = estBAdep - estBAdep(1,:).*factorq;
    estBAcum = estBAcum(2:end,:);
    estBAinst = IRLS_TV(estBAcum(:),P,muLocal,m-1,n,tol,[],ones((m-1)*n,1));
    estBAinst = reshape(estBAinst,m-1,n);

    %% ROI and metrics
    %  Masks
    xPm = medium.x; zPm = medium.z;
    [Xq, Zq] = meshgrid(xPm,zPm);
    Lx = radiusDisk*1; Lz = radiusDisk*1.3;
    cx = radiusDisk*1.45;

    inc = maskRect(xPm, zPm, 0, centerDepth, Lx, Lz);
    back = maskRect(xPm, zPm, -cx, centerDepth, Lx/2, Lz) | ...
        maskRect(xPm, zPm, cx, centerDepth, Lx/2, Lz);

    % Interp
    [Xmesh,Zmesh] = meshgrid(xP,zP(2:end));
    baInterp = interp2(Xmesh,Zmesh,estBAinst,Xq,Zq);
    acInterp = zeros(size(baInterp));

    % Metrics
    metricsDepRsld(iSim) = getMetrics(acInterp,baInterp,inc,back,'DepRsld', ...
        alphaInc/100,0,1);

    %% Plots

    figure('Position',imPosition);
    imagesc(xP*1e2,zP(2:end)*1e2,estBAinst); colorbar;
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

    %% ================================================================ %%
    %% Attenuation version 2
    % Plot and cumulative map
    attSam = (BR*filterParams.freqC)*100; % Np/cm to Np/m
    figure('Position',imPosition),
    imagesc(xP*100,zP*100,attSam)
    axis image
    c = colorbar;
    c.Label.String = 'dB/cm';
    title('Attenuation from RSLD')
    colormap turbo
    attSam = cumsum(attSam)./(1:length(zP))';
    estBAdep = (bz./( 1-exp(-2*attSam.*zP) ) .*attSam.*zP - 1)*2 ;

    %% Local maps with regularization
    [m,n] = size(bz);
    tol = 1e-3;
    muLocal = 0.01;
    dzP = zP(2)-zP(1);
    izP = round(zP./dzP);
    factorq = izP(1)./izP;
    P = sparse(tril(ones(m-1)));
    P = P./izP(2:end);
    P = kron(speye(n),P);

    estBAcum = estBAdep - estBAdep(1,:).*factorq;
    estBAcum = estBAcum(2:end,:);
    estBAinst = IRLS_TV(estBAcum(:),P,muLocal,m-1,n,tol,[],ones((m-1)*n,1));
    estBAinst = reshape(estBAinst,m-1,n);

    %% ROI and metrics
    %  Masks
    xPm = medium.x; zPm = medium.z;
    [Xq, Zq] = meshgrid(xPm,zPm);
    Lx = radiusDisk*1; Lz = radiusDisk*1.3;
    cx = radiusDisk*1.45;

    inc = maskRect(xPm, zPm, 0, centerDepth, Lx, Lz);
    back = maskRect(xPm, zPm, -cx, centerDepth, Lx/2, Lz) | ...
        maskRect(xPm, zPm, cx, centerDepth, Lx/2, Lz);

    % Interp
    [Xmesh,Zmesh] = meshgrid(xP,zP(2:end));
    baInterp = interp2(Xmesh,Zmesh,estBAinst,Xq,Zq);
    acInterp = zeros(size(baInterp));

    % Metrics
    metricsDepRsld2(iSim) = getMetrics(acInterp,baInterp,inc,back,'DepRsld2', ...
        alphaInc/100,0,1);

    %% Plots

    figure('Position',imPosition);
    imagesc(xP*1e2,zP(2:end)*1e2,estBAinst); colorbar;
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
T = struct2table([metricsDep,metricsDepRsld,metricsDepRsld2]);
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