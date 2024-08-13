%% Gauss-Newton
clear; close all; clc;
addpath(genpath(pwd))

nwl_vector = (5:5:50);
nwl_vector = 20;


%
% mu_vector = 10.^(-4:5);
% mu_vector = [10^-5 10^-1];
% %mu_vector = 10^-5;
% for idx =  1:length(mu_vector)
%    % disp(idx)
%
%     freq_vector = [5];
%
%     for ff = 1:length(freq_vector)
%         freq = freq_vector(ff);
%         mu = mu_vector(idx);
ff = 1; idx = 1;
freq = 5;
mu = 10;
nwl = nwl_vector;
% file_sam = ['rf_fnum3_PWNE_samBA12_att0p18f2_nc10_400kPa']
file_sam = ['rf_fnum3_PWNE_samBA9_att0p10f2_nc10_400kPa'];
file_ref = ['rf_fnum3_PWNE_refBA6_att10f2_nc10_400kPa'];

%alpha0_slope_dB = linspace(0.01,2,50)/freq; % dB/cm/MHz
%powerlaw = 2;
%      alpha0_slope_dB = linspace(0.0001,0.4,50); % dB/cm/MHz
%     powerlaw = 2;

filename = ['FULLMAPv60_',file_sam,'_',file_ref,'.mat'];

load(filename);
%%
BAmap = BAeff.image;
mzaxis = BAeff.mz;
z = BAeff.axial; % meters
x = BAeff.lateral;


dz = z(2)-z(1);
dx = x(2)-x(1);

% wl = 1540/1e6/freq;
wl = 1540/1e6/5;
blocksize = nwl*wl; % Change to test different sizes
nz = floor(blocksize/dz);
nx = floor(blocksize/dx);

zlim = [0 80]*1e-3;
zlim = [10 42]*1e-3;
%  zlim = [10 40]*1e-3;

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

while zini(end) < L1 + 1 - nz - param.nooverlap
    zini(end+1) = zini(end)+param.nooverlap;
end

zfin = zini + nz - 1;

nrow = length(zfin);
ncol = L2;



% alpha0_slope = alpha0_slope_dB/8.686; % Np/cm/MHz
f0 = freq;

%  alpha0 = alpha0_slope*f0^powerlaw*100;  % Np/m
ncol = L2;

%%
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
        betaR = 1+(10.5)/2;
        MLf = mzBA;%/betaR;

        ite = 1;
        % a = [a0 a1]
        while ite<100

            if ite == 1
                alphaArr(ite) = 0.1*5^2/8.686*100; % dB/cm -> Np/m
                betaArr(ite) = betaR;
            end

            dMLfda = -(betaArr(ite))*((2*zBA.*alphaArr(ite)+1).*exp(-2*alphaArr(ite)*zBA) - 1)./(alphaArr(ite)^2*zBA);
            dMLfdb = -(1./zBA/alphaArr(ite)).*(1 - exp(-2*alphaArr(ite)*zBA));

            jcb = [dMLfda dMLfdb];
            %mu = 0.01; %10^0;
            sqrtmu = sqrt(mu);
            jcb_lmm = [jcb; sqrtmu 0; 0 sqrtmu];


            res = ( MLf - (betaArr(ite)./zBA/alphaArr(ite)).*(1 - exp(-2*alphaArr(ite)*zBA)) );% residual

            res_lmm = [res; 0; 0];
            step = jcb\(-res);
            step = jcb_lmm\(-res_lmm);
            %keyboard
            alphaArr(ite+1) = alphaArr(ite) + step(1);
            betaArr(ite+1) = betaArr(ite) + step(2);
            %disp(['ite: ',ite,'. [alpha beta]: ',num2str(alpha(ite)/25/100*8.686), num2str(2*(beta(ite)-1)])

            norm_res(ite) = res'*res;
            ite = ite+1;
        end
        %figure; plot(norm_res);
        alphaArr;
        betaArr;
        %alpha_dB = alpha/5^1.05/100*8.686; % Np/cm/MH^2 : 0.1
        alpha_dB = alphaArr/5^2/100*8.686; % Np/cm/MH^2 : 0.1

        estAC(pixi,pixj) = alpha_dB(end);
        estBA(pixi,pixj) = 2*(betaArr(end)-1);
        clear alphaArr betaArr alpha_dB
    end
end

disp(['AC: ',num2str(nanmean(estAC(:))),' +/- ',num2str(nanstd(estAC(:)))]);
disp(['BA: ',num2str(nanmean(estBA(:))),' +/- ',num2str(nanstd(estBA(:)))]);
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
title(['B/A = ',num2str(nanmean(estBA(:)),'%.1f'),' +/- ',num2str(nanstd(estBA(:)),'%.1f')]);
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

mAC(idx) = mean2(estAC);
sAC(idx) = std2(estAC);
mBA(idx) = mean2(estBA);
sBA(idx) = std2(estBA);

%%

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
        betaR = 1+(10.5)/2;
        MLf = mzBA;%/betaR;

        ite = 1;
        alphaArr(ite) = 0.1*5^2/8.686*100; % dB/cm -> Np/m
        betaArr(ite) = betaR;

        while ite<100
            dMLfda = -(betaArr(ite))*((2*zBA.*alphaArr(ite)+1).*exp(-2*alphaArr(ite)*zBA) - 1)./(alphaArr(ite)^2*zBA);
            dMLfdb = -(1./zBA/alphaArr(ite)).*(1 - exp(-2*alphaArr(ite)*zBA));

            jcb = [dMLfda dMLfdb];
            res = ( MLf - (betaArr(ite)./zBA/alphaArr(ite)).*(1 - exp(-2*alphaArr(ite)*zBA)) );% residual
            step = jcb\(-res);

            alphaArr(ite+1) = alphaArr(ite) + step(1);
            betaArr(ite+1) = betaArr(ite) + step(2);
            %disp(['ite: ',ite,'. [alpha beta]: ',num2str(alpha(ite)/25/100*8.686), num2str(2*(beta(ite)-1)])

            norm_res(ite) = res'*res;
            ite = ite+1;
        end
        %figure; plot(norm_res);
        % alpha;
        % beta;
        %alpha_dB = alpha/5^1.05/100*8.686; % Np/cm/MH^2 : 0.1
        alpha_dB = alphaArr/5^2/100*8.686; % Np/cm/MH^2 : 0.1

        estAC(pixi,pixj) = alpha_dB(end);
        estBA(pixi,pixj) = 2*(betaArr(end)-1);
        % clear alphaArr betaArr alpha_dB
    end
end

disp(['AC: ',num2str(nanmean(estAC(:))),' +/- ',num2str(nanstd(estAC(:)))]);
disp(['BA: ',num2str(nanmean(estBA(:))),' +/- ',num2str(nanstd(estBA(:)))]);
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
title(['B/A = ',num2str(nanmean(estBA(:)),'%.1f'),' +/- ',num2str(nanstd(estBA(:)),'%.1f')]);
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

%%
