% Gauss-Newton
clear; close all; clc;


addpath(genpath(pwd))
nwl_vector = 40;
ff = 1; idx = 1;
freq = 5;
mu = 0;
nwl = nwl_vector;

%% Initialization
% file_sam = ['rf_fnum3_PWNE_samBA9_att0p10f2_nc10_400kPa'];
file_sam = 'rf_fnum3_PWNE_samBA9_att0p18f2_nc10_400kPa';
% file_sam = 'rf_fnum3_PWNE_samBA12_att0p18f2_nc10_400kPa';
file_ref = 'rf_fnum3_PWNE_refBA6_att10f2_nc10_400kPa';
beta0 = 1+(10.5)/2;
alpha0 = 0.18*5^2/8.686*100; % dB/cm -> Np/m

baseDir = 'C:\Users\smerino.C084288\Documents\MATLAB\BA_AC_joint\maps';
% baseDir = ['C:\Users\sebas\Documents\MATLAB\totalvarsimul_AC_BA\' ...
%      'BA_AC_joint\maps'];
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

%% Block by block
for pixj = 1:ncol
    for pixi = 1:nrow

        zblock = (z(zini(pixi)):1e-3:z(zfin(pixi))); % m
        clear index zBA mzBA
        for zi = 1:length(zblock)
            [~, idz] = min(abs(zblock(zi)-z));
            index(zi) = idz;
            zBA(zi,1) = z(idz);
            mzBA(zi,1) = mzaxis(idz,pixj);

        end
        MLf = mzBA;%/betaR;

        alpha(1) = alpha0;
        beta(1) = beta0;
        ite = 1;
        while ite<100
            dMLfda = -(beta(ite))*((2*zBA.*alpha(ite)+1).*exp(-2*alpha(ite)*zBA) - 1)./(alpha(ite)^2*zBA);
            dMLfdb = -(1./zBA/alpha(ite)).*(1 - exp(-2*alpha(ite)*zBA));

            jcb = [dMLfda dMLfdb];
            res = ( MLf - (beta(ite)./zBA/alpha(ite)).*(1 - exp(-2*alpha(ite)*zBA)) );% residual
            % step = jcb\(-res);
            [step,~] = cgs(jcb'*jcb,jcb'*-res);
            alpha(ite+1) = alpha(ite) + step(1);
            beta(ite+1) = beta(ite) + step(2);
            ite = ite+1;
        end
        alpha_dB = alpha/5^2/100*8.686; % Np/cm/MH^2 : 0.1

        estAC(pixi,pixj) = alpha_dB(end);
        estBA(pixi,pixj) = 2*(beta(end)-1);
        clear alpha beta alpha_dB
    end
end
disp(['AC: ',num2str(median(estAC(:), 'omitnan'))]);
disp(['BA: ',num2str(median(estBA(:), 'omitnan'))]);

%keyboard
figure; imagesc(x*1e3,z*1e3,estAC); colorbar;
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
figure; imagesc(x*1e3,z*1e3,estBA); colorbar;
clim([6 15]);
title('B/A');
title(['B/A = ',num2str(median(estBA(:),'omitnan'),'%.1f')]);
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

%% All at once (Coila's version)
zblock = (z(zini(1)):dz:z(zfin(1)));
p = length(zblock);

jcb_mtx = sparse(nrow*ncol*p, nrow*ncol*2);
res_mtx = sparse(nrow*ncol*p, 1);
jcb_lmm_mtx = sparse(nrow*ncol*(p+2), nrow*ncol*2);
res_lmm_mtx = sparse(nrow*ncol*(p+2), 1);
zBA_mtx = nan(nrow,ncol,p);
mzBA_mtx = nan(nrow,ncol,p);

for pixj = 1:ncol
    for pixi = 1:nrow

        zblock = (z(zini(pixi)):1e-3:z(zfin(pixi))); % m
        clear index zBA mzBA
        for zi = 1:length(zblock)
            zblock;
            % keyboard
            [~, idz] = min(abs(zblock(zi)-z));
            index(zi) = idz;
            zBA(zi,1) = z(idz);
            mzBA(zi,1) = mzaxis(idz,pixj);

        end
        zBA_mtx(pixi,pixj,:) = reshape(zBA, [1, 1, p]);
        mzBA_mtx(pixi,pixj,:) = reshape(mzBA, [1, 1, p]);
    end
end


ite = 1;

alpha_mtx(:) = alpha0; % dB/cm -> Np/m
beta_mtx(:) = beta0;
% permutIndex = ([1:ncol*nrow,1:ncol*nrow])';
% permutIndex(2:2:end) = permutIndex(2:2:end)+ncol*nrow;

while ite < 100
    for pixj = 1:ncol
        for pixi = 1:nrow

            zBA = reshape( zBA_mtx(pixi,pixj,:) , [p, 1 ,1]);
            mzBA = reshape( mzBA_mtx(pixi,pixj,:) , [p, 1 ,1]);
            MLf = mzBA;%/beta0;

            beta(ite) = beta_mtx(nrow,ncol);
            alpha(ite) = alpha_mtx(nrow,ncol);

            dMLfda = -(beta(ite))*((2*zBA.*alpha(ite)+1).*exp(-2*alpha(ite)*zBA) - 1)./(alpha(ite)^2*zBA);
            dMLfdb = -(1./zBA/alpha(ite)).*(1 - exp(-2*alpha(ite)*zBA));

            jcb = [dMLfda dMLfdb];
            res = ( MLf - (beta(ite)./zBA/alpha(ite)).*(1 - exp(-2*alpha(ite)*zBA)) );% residual

            idx = pixi + (pixj-1)*nrow; % it need to go from 1 to nrow*ncol
            jcb_mtx( idx*p - (p-1): idx*p , idx*2 - 1: idx*2 ) = jcb;
            res_mtx( idx*p - (p-1): idx*p , 1 ) = res;

            % Levenberg-Maquard version
            % sqrtmu = sqrt(mu);

            % jcb_lmm = [jcb; sqrtmu 0; 0 sqrtmu];
            % res_lmm = [res; 0; 0];
            % jcb_lmm_mtx( idx*(p+2) - (p+2-1): idx*(p+2) , idx*2 - 1: idx*2 ) = jcb_lmm;
            % res_lmm_mtx( idx*(p+2) - (p+2-1): idx*(p+2) , 1 ) = res_lmm;

        end
    end
    % jcb_mtx = jcb_mtx(:,[1:2:ncol*nrow*2, 2:2:ncol*nrow*2]);

    % ab_mtx = cgs(jcb_lmm_mtx'*jcb_lmm_mtx,jcb_lmm_mtx'*(-res_lmm_mtx));
    [ab_mtx,~] = cgs(jcb_mtx'*jcb_mtx,jcb_mtx'*(-res_mtx));
    % ab_mtx = ab_mtx(permutIndex);

    % ab_mtx = pcg(jcb_mtx'*jcb_mtx,jcb_mtx'*(-res_mtx));
    step_alpha_mtx_int = ab_mtx(1:2:end);
    step_alpha_mtx = reshape(step_alpha_mtx_int, [nrow, ncol]);
    step_beta_mtx_int = ab_mtx(2:2:end);
    step_beta_mtx = reshape(step_beta_mtx_int, [nrow, ncol]);

    alpha_mtx = alpha_mtx + step_alpha_mtx;
    beta_mtx = beta_mtx + step_beta_mtx;

    norm_res(ite) = res_mtx'*res_mtx;
    ite = ite+1;
end

alpha_dB_mtx = alpha_mtx/5^2/100*8.686; % Np/cm/MH^2 : 0.1 or 0.18

estACfail = alpha_dB_mtx;
estBAfail = 2*(beta_mtx-1);

disp(['AC: ',num2str(median(estACfail(:), 'omitnan'))]);
disp(['BA: ',num2str(median(estBAfail(:), 'omitnan'))]);

figure; imagesc(x*1e3,z*1e3,estACfail); colorbar;
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
figure; imagesc(x*1e3,z*1e3,estBAfail); colorbar;
clim([6 15]);
title('B/A');
title(['B/A = ',num2str(median(estBAfail(:),'omitnan'),'%.1f')]);
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
for ite =1:100
    jcb = jacobian(theta,X);
    res = Y - model(theta,X);
    [step,~] = cgs(jcb'*jcb + mu*speye(2*ncol*nrow),jcb'*-res);
    % step = (jcb'*jcb)\(jcb'*-res);
    theta = theta + step;
end

alphaArr = theta(1:ncol*nrow);
betaArr = theta(ncol*nrow+1:end);
alpha_dB = alphaArr/5^2/100*8.686; % Np/cm/MH^2 : 0.1

estAC2 = reshape(alpha_dB,[nrow,ncol]);
estBA2 = reshape(2*(betaArr-1),[nrow,ncol]);

disp(['AC: ',num2str(median(estAC2(:), 'omitnan'))]);
disp(['BA: ',num2str(median(estBA2(:), 'omitnan'))]);

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
title(['B/A = ',num2str(median(estBA2(:),'omitnan'),'%.1f')]);
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

% save_all_figures_to_directory('figures','fig','fig')

fprintf('\nMSE: %.3f\n\n',mean((estBA(:)-estBA2(:)).^2))

%% My version two
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
theta = repmat([alpha0;beta0],1,ncol*nrow);
for ite =1:100
    jcb = jacobian(theta,X);
    res = Y - model(theta,X);
    [step,~] = cgs(jcb'*jcb + mu*speye(2*ncol*nrow),jcb'*-res);
    % step = (jcb'*jcb)\(jcb'*-res);
    theta = theta + step;
end

alphaArr = theta(1:ncol*nrow);
betaArr = theta(ncol*nrow+1:end);
alpha_dB = alphaArr/5^2/100*8.686; % Np/cm/MH^2 : 0.1

estAC2 = reshape(alpha_dB,[nrow,ncol]);
estBA2 = reshape(2*(betaArr-1),[nrow,ncol]);

disp(['AC: ',num2str(median(estAC2(:), 'omitnan'))]);
disp(['BA: ',num2str(median(estBA2(:), 'omitnan'))]);

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
title(['B/A = ',num2str(median(estBA2(:),'omitnan'),'%.1f')]);
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


%% Comparing first iteration
% Initialization
alpha_mtx(:) = alpha0; % dB/cm -> Np/m
beta_mtx(:) = beta0;
ite = 1;
theta = [alpha0*ones(ncol*nrow,1);beta0*ones(ncol*nrow,1)]; % Mine


% Coila's version
for pixj = 1:ncol
    for pixi = 1:nrow

        zBA = reshape( zBA_mtx(pixi,pixj,:) , [p, 1 ,1]);
        mzBA = reshape( mzBA_mtx(pixi,pixj,:) , [p, 1 ,1]);
        MLf = mzBA;%/beta0;

        beta(ite) = beta_mtx(nrow,ncol);
        alpha(ite) = alpha_mtx(nrow,ncol);

        dMLfda = -(beta(ite))*((2*zBA.*alpha(ite)+1).*exp(-2*alpha(ite)*zBA) - 1)./(alpha(ite)^2*zBA);
        dMLfdb = -(1./zBA/alpha(ite)).*(1 - exp(-2*alpha(ite)*zBA));

        jcb = [dMLfda dMLfdb];
        res = ( MLf - (beta(ite)./zBA/alpha(ite)).*(1 - exp(-2*alpha(ite)*zBA)) );% residual

        idx = pixi + (pixj-1)*nrow; % it need to go from 1 to nrow*ncol
        jcb_mtx( idx*p - (p-1): idx*p , idx*2 - 1: idx*2 ) = jcb;
        res_mtx( idx*p - (p-1): idx*p , 1 ) = res;
    end
end
[ab_mtx,~] = cgs(jcb_mtx'*jcb_mtx,jcb_mtx'*(-res_mtx));
step_alpha_mtx_int = ab_mtx(1:2:end);
step_alpha_mtx = reshape(step_alpha_mtx_int, [nrow, ncol]);
step_beta_mtx_int = ab_mtx(2:2:end);
step_beta_mtx = reshape(step_beta_mtx_int, [nrow, ncol]);

alpha_mtx = alpha_mtx + step_alpha_mtx;
beta_mtx = beta_mtx + step_beta_mtx;

norm_res(ite) = res_mtx'*res_mtx;
ite = ite+1;


% My version
jcb = jacobian(theta,X);
res = Y - model(theta,X);
[step,~] = cgs(jcb'*jcb + mu*speye(2*ncol*nrow),jcb'*-res);
theta = theta + step;


jcb_mtx_permuted = jcb_mtx(:,[1:2:ncol*nrow*2, 2:2:ncol*nrow*2]);
ab_mtx_permuted = ab_mtx([1:2:ncol*nrow*2, 2:2:ncol*nrow*2]);
disp('Are Jacobians equal?')
disp(all(jcb(:) == jcb_mtx_permuted(:)))
disp('Are residuals equal?')
disp(all(res(:) == res_mtx(:)))
disp('Are steps equal?')
disp(all(step(:) == ab_mtx_permuted(:)))
disp(norm(ab_mtx_permuted - step))
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


function jcb = jacobian2(theta,x)
mn = size(theta,1)/2; % number of points in image
p = length(x)/mn;

alphaArr = theta(1:2:end)';
betaArr = theta(1:2:end)';
zBA = reshape(x,[p,mn]);

dMLfda = -(betaArr).*((2*zBA.*alphaArr+1).*exp(-2*alphaArr.*zBA) - 1)./(alphaArr.^2.*zBA);
dMLfdb = -(1./zBA./alphaArr).*(1 - exp(-2*alphaArr.*zBA));


% Indices for sparse matrix
I1 = reshape(1:2:2*mn*p,[mn,p]); % Row indices
I2 = reshape(2:2:2*mn*p,[mn,p]); % Row indices

J = (1:mn).*ones(p,1);

jcb = [sparse([I1;I2],[J;J],[dMLfda;dMLfdb])];
end

function mdl = model2(theta,x)
mn = size(theta,1)/2; % number of points in image
p = length(x)/mn;

alphaArr = theta(1:2:end)';
betaArr = theta(1:2:end)';
zBA = reshape(x,[p,mn]);

mdl = (betaArr./zBA./alphaArr).*(1 - exp(-2*alphaArr.*zBA)) ;
mdl = mdl(:);
end
