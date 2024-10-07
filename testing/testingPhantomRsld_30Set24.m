setup

dataDir = 'C:\Users\sebas\Documents\Data\Nonlinearity\30-Set-24\bf_wideBand';
refDir = 'C:\Users\sebas\Documents\Data\Nonlinearity\adriana2';
% refDir = dataDir;
resultsDir = fullfile(dataDir,'results_2');
[~,~,~] = mkdir(resultsDir);
targetFiles = dir([dataDir,'\*5MHz*_2.mat']);

%%
blocksize = 16;     % Block size in wavelengths
freq_L = 4e6; freq_H = 7e6;
freq_C = (freq_L + freq_H)/2;

overlap_pc      = 0.8;
ratio_zx        = 3/2;
NptodB = log10(exp(1))*20;

% Weight parameters
muB = 10^3; muC = 10^0;
ratioCutOff = 10;
order = 5;
reject = 0.1;
extension = 3;

% SWTV
aSNR = 5; bSNR = 0.09;
desvMin = 15;

groundTruthTargets = [0.97,0.95,0.95,0.55];

% Plotting constants
dynRange = [-50,0];
attRange = [0,2.5];

tol = 1e-3;

c1x = 1.95; c1z = 1.93;
roiL = 1; roiD = 0.6;
roiLz = 1.5;

% muBtv = 10^3.5; muCtv = 10^3.5;
muBwfr = 10^3.5; muCwfr = 10^2;
muBtv = 10^2.5; muCtv = 10^2.5;

iAcq = 3;

%%
for iAcq = 1:length(targetFiles)
fprintf("Simulation no. %i, %s\n",iAcq,targetFiles(iAcq).name);
out = load(fullfile(dataDir,targetFiles(iAcq).name));
dx = out.x(2)-out.x(1);
dz = out.z(2)-out.z(1);
x = out.x*1e2; % [cm]
z = out.z'*1e2; % [cm]
c0 = out.c0;
sam1 = out.rf1;
fs = out.fs;

%% Cropping and finding sample sizes
x_inf = -2; x_sup = 2;
z_inf = 0.2; z_sup = 3;

% Limits for ACS estimation
ind_x = x_inf <= x & x <= x_sup;
ind_z = z_inf <= z & z <= z_sup;
roi = ind_x.*ind_z;
x = x(ind_x);
z = z(ind_z);
sam1 = sam1(ind_z,ind_x);

% Wavelength size
wl = c0/mean(freq_C);   % Wavelength (m)

% Lateral samples
wx = round(blocksize*wl*(1-overlap_pc)/dx);  % Between windows
nx = round(blocksize*wl/dx);                 % Window size
x0 = 1:wx:length(x)-nx;
x_ACS = x(1,x0+round(nx/2));
n  = length(x0);

% Axial samples
wz = round(blocksize*wl*(1-overlap_pc)/dz * ratio_zx); % Between windows
nz = 2*round(blocksize*wl/dz /2 * ratio_zx); % Window size
L = (nz/2)*dz*100;   % (cm)
z0p = 1:wz:length(z)-nz;
z0d = z0p + nz/2;
z_ACS = z(z0p+ nz/2);
m  = length(z0p);

[pxx,fpxx] = pwelch(sam1-mean(sam1),500,400,500,fs);
meanSpectrum = db(mean(pxx,2));
meanSpectrum = meanSpectrum - max(meanSpectrum);
figure,plot(fpxx/1e6,meanSpectrum)
% [freq_L,freq_H] = findFreqBand(fpxx, meanSpectrum, 0.1);
xline([freq_L,freq_H]/1e6)
xlim([0 fs/2/1e6])
xlabel('Frequency [MHz]')
ylabel('Magnitude [dB]')
grid on
title('Frequency spectrum')

% Frequency samples
NFFT = 2^(nextpow2(nz/2)+2);
band = (0:NFFT-1)'/NFFT * fs;   % [Hz] Band of frequencies
rang = band > freq_L & band < freq_H ;   % useful frequency range
f  = band(rang)*1e-6; % [MHz]
p = length(f);

% Plot region of interest B-mode image
Bmode = db(hilbert(sam1));
Bmode = Bmode - max(Bmode(:));

fprintf('\nFrequency range: %.2f - %.2f MHz\n',freq_L*1e-6,freq_H*1e-6)
fprintf('Blocksize in wavelengths: %i\n',blocksize)
fprintf('Blocksize x: %.2f mm, z: %.2f mm\n',nx*dx*1e3,nz*dz*1e3)
fprintf('Blocksize in pixels nx: %i, nz: %i\n',nx,nz);
fprintf('Region of interest columns: %i, rows: %i\n\n',m,n);



%% Generating Diffraction compensation
% Generating references
att_ref = 0.5*f.^1.12/NptodB; 
att_ref_map = zeros(m,n,p);
for jj=1:n
    for ii=1:m
        att_ref_map(ii,jj,:) = att_ref;
    end
end

% Windows for spectrum
windowing = tukeywin(nz/2,0.25);
% windowing = hamming(nz/2);
windowing = windowing*ones(1,nx);

% For looping
refFiles = dir(fullfile(refDir,'*.mat'));
Nref = length(refFiles);

% Memory allocation
Sp_ref = zeros(m,n,p,Nref);
Sd_ref = zeros(m,n,p,Nref);
for iRef = 1:Nref
    out = load(fullfile(refDir,refFiles(iRef).name));
    samRef = out.ps_datal(:,:,1);
    samRef = samRef(ind_z,ind_x); % Cropping
    % figure,imagesc(db(hilbert(samRef)))
    for jj=1:n
        for ii=1:m
            xw = x0(jj) ;   % x window
            zp = z0p(ii);
            zd = z0d(ii);

            sub_block_p = samRef(zp:zp+nz/2-1,xw:xw+nx-1);
            sub_block_d = samRef(zd:zd+nz/2-1,xw:xw+nx-1);
            [tempSp,~] = spectra(sub_block_p,windowing,0,nz/2,NFFT);
            [tempSd,~] = spectra(sub_block_d,windowing,0,nz/2,NFFT);

            Sp_ref(ii,jj,:,iRef) = (tempSp(rang));
            Sd_ref(ii,jj,:,iRef) = (tempSd(rang));
        end
    end
end

Sp = mean(Sp_ref,4); Sd = mean(Sd_ref,4);
compensation = ( log(Sp) - log(Sd) ) - 4*L*att_ref_map;

% Liberating memory to avoid killing my RAM
clear Sp_ref Sd_ref

%% Spectrum
Sp = zeros(m,n,p);
Sd = zeros(m,n,p);
for jj=1:n
    for ii=1:m
        xw = x0(jj) ;   % x window
        zp = z0p(ii);
        zd = z0d(ii);

        sub_block_p = sam1(zp:zp+nz/2-1,xw:xw+nx-1);
        sub_block_d = sam1(zd:zd+nz/2-1,xw:xw+nx-1);

        [tempSp,~] = spectra(sub_block_p,windowing,0,nz/2,NFFT);
        [tempSd,~] = spectra(sub_block_d,windowing,0,nz/2,NFFT);
        Sp(ii,jj,:) = (tempSp(rang));
        Sd(ii,jj,:) = (tempSd(rang));
    end
end

%% ROI selection
[X,Z] = meshgrid(x_ACS,z_ACS);
[Xq,Zq] = meshgrid(x,z);
% rInc = 0.95;
% inc = ((Xq-c1x).^2 + (Zq-c1z).^2)<= (rInc-0.1)^2;
% back = ((Xq-c1x).^2 + (Zq-c1z).^2) >= (rInc+0.1)^2;
% 
% x0mask = c1x - roiL/2; 
% z0mask = c1z - roiLz/2;
% [back,inc] = getRegionMasks(x,z,c1x,c1z,roiL,roiD,roiLz);

% Setting up
b = (log(Sp) - log(Sd)) - (compensation);
A1 = kron( 4*L*f , speye(m*n) );
A2 = kron( ones(size(f)) , speye(m*n) );
mask = ones(m,n,p);

%% RSLD-TV

tic
[Bn,Cn] = AlterOpti_ADMM(A1,A2,b(:),muBtv,muCtv,m,n,tol,mask(:));
toc
BR = (reshape(Bn*NptodB,m,n));
CR = (reshape(Cn*NptodB,m,n));

%% SWIFT
[~,Cn] = optimAdmmTvTikhonov(A1,A2,b(:),muBwfr,muCwfr,m,n,tol,mask(:));
bscMap = reshape(Cn,m,n)*NptodB;

w = (1-reject)* (1./((bscMap/ratioCutOff).^(2*order) + 1)) + reject;
wExt = movmin(w,extension);

% New equations
W = repmat(wExt,[1 1 p]);
W = spdiags(W(:),0,m*n*p,m*n*p);
bw = W*b(:);
A1w = W*A1;
A2w = W*A2;

% Regularization: Au = b
tic
[Bn,Cn] = optimAdmmWeightedTvTikhonov(A1w,A2w,bw,muBwfr,muCwfr,m,n,tol,mask(:),w);
toc
BSWIFT = (reshape(Bn*NptodB,m,n));
CSWIFT = (reshape(Cn*NptodB,m,n));


%% Figures
c1x = (-1.45+1.6)/2 + 0.1;
c1z = (2.95+5.8)/2;
rInc = (5.8-2.95)/4 + (1.45+1.6)/4;

[X,Z] = meshgrid(x_ACS,z_ACS);
% regionMaskAcs = (X - c1x).^2 + (Z-c1z).^2 < rInc^2;
regionMaskAcs = true(m,n);

acsTv = sum(BR.*regionMaskAcs,"all")/sum(regionMaskAcs,"all");
zcTv = sum(CR.*regionMaskAcs,"all")/sum(regionMaskAcs,"all")/4/L;
acsSwift = sum(BSWIFT.*regionMaskAcs,"all")/sum(regionMaskAcs,"all");
zcSwift = sum(CSWIFT.*regionMaskAcs,"all")/sum(regionMaskAcs,"all")/4/L;


figure('Units','centimeters', 'Position',[5 5 14 7]);
tl = tiledlayout(1,2, "Padding","tight");

t1 = nexttile;
imagesc(x,z,Bmode,dynRange)
xlim([x_ACS(1) x_ACS(end)]),
ylim([z_ACS(1) z_ACS(end)]),
xlabel('Lateral [cm]'),
ylabel('Axial [cm]')
axis image
colormap(t1,gray)
title('B-mode')
%subtitle(' ')
c = colorbar;
c.Label.String = 'dB';
% hold on 
% rectangle('Position',[c1x-rInc c1z-rInc 2*rInc 2*rInc], 'LineStyle','--', ...
%     'LineWidth',1, 'Curvature',1)
% hold off

% fontsize(gcf,8,'points')

t2 = nexttile;
imagesc(x_ACS,z_ACS,BR, attRange)
xlabel('Lateral [cm]'),
ylabel('Axial [cm]')
colormap(t2,turbo)
axis image
title(sprintf('RSLD, ACS = %.2f',acsTv))
c = colorbar;
c.Label.String = 'dB/cm/MHz';
% fontsize(gcf,8,'points')
% hold on 
% rectangle('Position',[c1x-rInc c1z-rInc 2*rInc 2*rInc], 'LineStyle','--', ...
%     'LineWidth',1, 'Curvature',1)
% hold off


% t5 = nexttile;
% imagesc(x_ACS,z_ACS,BSWIFT, attRange)
% xlabel('Lateral [cm]'),
% ylabel('Axial [cm]')
% colormap(t5,turbo)
% axis image
% title(sprintf('SWIFT, ACS = %.2f',acsSwift))
% c = colorbar;
% c.Label.String = 'dB/cm/MHz';

% fontsize(gcf,8,'points')
% hold on 
% rectangle('Position',[c1x-rInc c1z-rInc 2*rInc 2*rInc], 'LineStyle','--', ...
%     'LineWidth',1, 'Curvature',1)
% hold off
%% SLD fit and Weighted SLD fit

sldLine = squeeze(sum(sum(b.*regionMaskAcs,2),1))/4/L*NptodB/sum(regionMaskAcs(:));

figure('Units','centimeters', 'Position',[5 5 10 10]),
% tiledlayout(1,2),
% nexttile
plot(f,sldLine),
hold on,
plot(f,acsTv*f + zcTv, '--')
hold off,
grid on,
axis([0 freq_H/1e6 0 16])
xlabel('Frequency [MHz]')
ylabel('Attenuation [dB/cm]')
title('Mean SLD')
legend('SLD',sprintf('AC=%.2f%+.2ff\n',zcTv,acsTv), 'Location','northwest')

%% Power law fit
fit2 = [ones(length(f),1) log(f)]\log(sldLine);
alphaCoeff = exp(fit2(1)); alphaPower = fit2(2);
figure('Units','centimeters', 'Position',[5 5 10 10]),
% tiledlayout(1,2),
% nexttile
plot(f,sldLine),
hold on,
plot(f,alphaCoeff*f.^alphaPower,'--')
hold off,
grid on,
axis([0 freq_H/1e6 0 16])
xlabel('Frequency [MHz]')
ylabel('Attenuation [dB/cm]')
title('Power law fit')
legend('SLD',sprintf('AC=%.2ff^{%.2f}\n',alphaCoeff,alphaPower), 'Location','northwest')

%% Square fit
fit2 = (f.^2)\sldLine;
alphaCoeff = (fit2(1));
figure('Units','centimeters', 'Position',[5 5 10 10]),
% tiledlayout(1,2),
% nexttile
plot(f,sldLine),
hold on,
plot(f,alphaCoeff*f.^2,'--')
hold off,
grid on,
axis([0 freq_H/1e6 0 16])
xlabel('Frequency [MHz]')
ylabel('Attenuation [dB/cm]')
title('Power law fit')
legend('SLD',sprintf('AC=%.3ff^{2}\n',alphaCoeff), 'Location','northwest')


%%
save_all_figures_to_directory(resultsDir,char("sam"+iAcq+"fig"),'fig')
close all
end
