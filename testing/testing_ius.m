% Gauss-Newton
clear; close all; clc;
addpath(genpath(pwd))
nwl_vector = 40;
baseDir = 'C:\Users\smerino.C084288\Documents\MATLAB\BA_AC_joint\maps';
%%
mu_vector = 10.^(-4:1:5);
%mu_vector = [10^-5 10^-1];
%mu_vector = 10^-5;
for iMu =  1:length(mu_vector)
    %% Hyperparameters
    file_sam = 'rf_fnum3_PWNE_samBA9_att0p10f2_nc10_400kPa';
    % file_sam = 'rf_fnum3_PWNE_samBA12_att0p18f2_nc10_400kPa';

    file_ref = 'rf_fnum3_PWNE_refBA6_att10f2_nc10_400kPa';

    beta0 = 1+(10.5)/2;
    alpha0 = 0.08*5^2/8.686*100; % dB/cm -> Np/m
    resultName = fullfile('results','muBdivmuA1','samBA9_attIni0p08');

    freq = 5;
    muBeta = mu_vector(iMu); muAlpha = muBeta/1;
    nwl = nwl_vector;

    %% Loading and cropping
    filename = ['FULLMAPv60_',file_sam,'_',file_ref,'.mat'];

    load(fullfile(baseDir,filename));
    BAmap = BAeff.image;
    mzaxis = BAeff.mz;
    z = BAeff.axial; % meters
    x = BAeff.lateral;


    dz = z(2)-z(1);
    dx = x(2)-x(1);

    wl = 1540/1e6/freq;
    blocksize = nwl*wl; % Change to test different sizes
    nz = floor(blocksize/dz);
    nx = floor(blocksize/dx);

    % zlim = [0 80]*1e-3;
    zlim = [10 42]*1e-3;

    zori = z;
    mzaxisori = mzaxis;
    [~,idzmin] = min(abs(zlim(1)-z));
    [~,idzmax] = min(abs(zlim(2)-z));
    z = zori(idzmin:idzmax);
    mzaxis = mzaxisori(idzmin:idzmax,:);
    L1   = size(mzaxis,1);
    L2   = size(mzaxis,2);

    param.overlap_pct = 0.1;
    param.nz = nz;
    param.overlap = round((1-param.overlap_pct)*param.nz);
    param.nooverlap = nz - param.overlap;
    zini(1)  = 1;
    dz = 1e-3;

    while zini(end) < L1 + 1 - nz - param.nooverlap
        zini(end+1) = zini(end)+param.nooverlap;
    end

    zfin = zini + nz - 1;
    f0 = freq;
    nrow = length(zfin);
    ncol = L2;

    %% MY VERSION
    % Constructing arrays
    zblock = (z(zini(1)):1e-3:z(zfin(1))); % m
    np = length(zblock); % Number of points per block

    mzBA = zeros(np,nrow*ncol);
    zBA = zeros(np,nrow*ncol); % zBA
    for pixj = 1:ncol
        for pixi = 1:nrow
            zblock = (z(zini(pixi)):1e-3:z(zfin(pixi))); % m
            for zi = 1:length(zblock)
                [~, idz] = min(abs(zblock(zi)-z));
                zBA(zi,pixi+(pixj-1)*nrow) = z(idz);
                mzBA(zi,pixi+(pixj-1)*nrow) = mzaxis(idz,pixj);
                % X((pixj-1)*np+(pixi-1)*ncol*np+zi) = z(idz);
                % Y((pixj-1)*np+(pixi-1)*ncol*np+zi) = mzaxis(idz,pixj);
            end
        end
    end

    X = zBA(:); % zBA
    Y = mzBA(:); % mzBA

    % Iterating
    theta = [alpha0*ones(ncol*nrow,1);beta0*ones(ncol*nrow,1)];
    regMatrix = blkdiag(muAlpha*speye(ncol*nrow),muBeta*speye(ncol*nrow));
    for ite =1:100
        jcb = jacobian(theta,X);
        res = Y - model(theta,X);
        [step,~] = pcg(jcb'*jcb + regMatrix,jcb'*-res);
        % step = (jcb'*jcb)\(jcb'*-res);
        theta = theta + step;
    end

    alphaArr = theta(1:ncol*nrow);
    betaArr = theta(ncol*nrow+1:end);
    alpha_dB = alphaArr/freq^2/100*8.686; % Np/cm/MH^2 : 0.1

    estAC2 = reshape(alpha_dB,[nrow,ncol]);
    estBA2 = reshape(2*(betaArr-1),[nrow,ncol]);

    disp(['AC: ',num2str(median(estAC2(:), 'omitnan'))]);
    disp(['BA: ',num2str(median(estBA2(:), 'omitnan'))]);
    
    %%
    figure; imagesc(x*1e3,z*1e3,estAC2); colorbar;
    clim([0 0.2]);
    %title('ACS (dB/cm/MHz)');
    title('\alpha_0 in \alpha(f) = \alpha_0 \times f^2 dB/cm')
    font = 20;
    axis image
    colormap pink; colorbar;
    %clim([0 12])
    %font = 20;
    set(gca,'fontsize',font)
    xlabel('Lateral distance (mm)');
    ylabel('Depth (mm)');
    %title("\bf BoA map");
    set(gca,'FontSize',font);
    %ylim([5 55])
    pause(0.5)
    figure; imagesc(x*1e3,z*1e3,estBA2); colorbar;
    clim([6 15]);
    title('B/A');
    title(['B/A = ',num2str(median(estBA2(:),'omitnan'),'%.1f'),...
        ', log_{10}\mu_{\beta} = ', num2str(log10(muBeta),'%i')]);
    font = 20;
    axis image
    colormap pink; colorbar;
    %clim([0 12])
    %font = 20;
    set(gca,'fontsize',font)
    xlabel('Lateral distance (mm)');
    ylabel('Depth (mm)');
    %title("\bf BoA map");
    set(gca,'FontSize',font);
    alpha_mtx = nan(nrow,ncol);
    beta_mtx = nan(nrow,ncol);

    %%
    % save_all_figures_to_directory('figures','LMfig','fig')
    % 
    % disp('MSE')
    % disp(mean((estBA(:)-estBA2(:)).^2))



    %keyboard
    mAC(iMu) = mean(estAC2(:),'omitnan');
    sAC(iMu) = std(estAC2(:),[],'omitnan');
    mBA(iMu) = mean(estBA2(:),'omitnan');
    sBA(iMu) = std(estBA2(:),[],'omitnan');

end

%%
figure
errorbar(mu_vector,mBA, sBA, 'r', 'LineWidth',3)
set(gca, 'XScale', 'log')
hold on;
plot(mu_vector,9*ones(size(mu_vector)),'k--','LineWidth',3)
%ylim([2 20])
%xlim([10^(-3.5) 10^(2.5)])
xlabel('\mu ')
title('B/A')
set(gca,'FontSize',20)
% 
save(resultName, 'mBA', 'sBA', 'mu_vector')

%% Utility functions
function jcb = jacobian(theta,x)
mn = size(theta,1)/2; % number of points in image
p = length(x)/mn;

alphaArr = theta(1:mn)';
betaArr = theta(mn+1:end)';
zBA = reshape(x,[p,mn]);

dMLfda = -(betaArr).*((2*zBA.*alphaArr+1).*exp(-2*alphaArr.*zBA) - 1)./(alphaArr.^2.*zBA);
dMLfdb = -(1./zBA./alphaArr).*(1 - exp(-2*alphaArr.*zBA));


% Indices for sparse matrix
I = reshape(1:mn*p,[mn,p]); % Row indices
J = (1:mn).*ones(p,1);

jcb = [sparse(I,J,dMLfda) sparse(I,J,dMLfdb)];
end

function mdl = model(theta,x)
mn = size(theta,1)/2; % number of points in image
p = length(x)/mn;

alphaArr = theta(1:mn)';
betaArr = theta(mn+1:end)';
zBA = reshape(x,[p,mn]);

mdl = (betaArr./zBA./alphaArr).*(1 - exp(-2*alphaArr.*zBA)) ;
mdl = mdl(:);
end