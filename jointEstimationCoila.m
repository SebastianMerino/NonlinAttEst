%% B/A estimator with time domain approach
% by: Andres Coila

clear; close all; clc;
addpath(genpath(pwd))


%% Prepare files for sanmple and reference

param.instfreq = 'no';   % 'yes' or 'no'


samDir = 'C:\Users\sebas\Documents\Data\Nonlinearity\attInc\fn2';
fileSam = "RFfn2_PWNE5MHz_samincBA6inc12_att0p1f2inc0p10_nc10_400kPa";
refDir = 'C:\Users\sebas\Documents\Data\Nonlinearity\uniformBA_inclusionACS';
fileRef = "RFfn2_PWNE5MHz_samBA12_att0p1f2inc0p10_nc10_400kPa";


load(fullfile(samDir,fileSam))
rf = rf1; save('rf_sam_LP','rf','fs','c0','x','z')
rf = rf2; save('rf_sam_HP','rf','fs','c0','x','z')
clear rf1 rf2
load(fullfile(refDir,fileRef))
rf = rf1; save('rf_ref_LP','rf','fs','c0','x','z')
rf = rf2; save('rf_ref_HP','rf','fs','c0','x','z')

media.samLP = 'rf_sam_LP';
media.samHP = 'rf_sam_HP';
media.refLP = 'rf_ref_LP';
media.refHP = 'rf_ref_HP';


%% Assign attenuation corresponding to the RF data
NptodB = 20*log10(exp(1));

param.ACsam = 0.1/NptodB;
param.APsam = 2;
param.ACref = 0.1/NptodB;
param.APref = 2;
param.BoAref = 12;

% [param.ACsam, param.APsam] = acs_filename(file_sam);
% [param.ACref, param.APref] = acs_filename(file_ref);
% param.BoAref = BonA_filename(file_ref);

%% Asign other hyperparameters. Double check!

param.factor = 5;   % when P goes from 100 kPa to 1000 kPa
param.order_filter = 200; % Filter for time domain rf data
% param.order_filter = 50; % Filter for time domain rf data
param.f0 = 5;   % MHz
param.size_wl = 10;%4.5; %4;%; %10;   % Moving average window for envelopes

param.width = 4; % Width to smooth envelopes laterally per B/A line
param.overlap_pct = 0.5;    % 80%
param.overlap = round((1-param.overlap_pct)*param.width);

param.col_last = param.width;
param.n = 0;

param.fs = fs;
param.fnumber = 3; % if plane wave Beamformingstill required

param.noise_level = 0;
%param.alpha_power = 1.05;

%% Estimation of B/A image

[BAimage, BAparam] = BAestimator(media, param);

%%
figure; imagesc(BAparam.x, BAparam.z, BAimage);
font = 20;
axis image
colormap pink; colorbar;
clim([5 10])
%font = 20;
set(gca,'fontsize',font)
xlabel('\bfLateral distance (mm)');
ylabel('\bfDepth (mm)');
%title("\bf BoA map");
set(gca,'FontSize',font);
ylim([5 55])

figure; plot(BAparam.z, BAimage(:,ceil(size(BAimage,2)/2)),'LineWidth',3);
grid on;
set(gca,'fontsize',font)
xlabel('\bfDepth (mm)');
ylabel('\bf B/A');
%  title("\bf B/A");
set(gca,'FontSize',font);
ylim([0 15])


%keyboard
%return

BAeff.image = BAimage; %QUS
BAeff.lateral = BAparam.x*1e-3;
BAeff.axial = BAparam.z*1e-3;
BAeff.mz = BAparam.mz;


%% Getting system
freq = 5;
muAlpha = 1; muBeta = 1;
zlim = [5 55]*1e-3;
nwl = 40;
beta0 = 1+(10.5)/2;
alpha0 = 0.1*5^2/8.686*100; % dB/cm -> Np/m

% mzaxis = movmean(BAeff.mz,4,2);
mzaxis = BAeff.mz;
x = BAeff.lateral';
dz = z(2)-z(1);
dx = x(2)-x(1);

wl = 1540/1e6/freq;
blocksize = nwl*wl; % Change to test different sizes
nz = floor(blocksize/dz);
nx = floor(blocksize/dx);


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
m = length(zfin);
n = L2;

%% MY VERSION
% Constructing arrays
zblock = (z(zini(1)):1e-3:z(zfin(1))); % m
np = length(zblock); % Number of points per block

mzBA = zeros(np,m*n);
zBA = zeros(np,m*n); % zBA
for pixj = 1:n
    for pixi = 1:m
        zblock = (z(zini(pixi)):1e-3:z(zfin(pixi))); % m
        for zi = 1:length(zblock)
            [~, idz] = min(abs(zblock(zi)-z));
            zBA(zi,pixi+(pixj-1)*m) = z(idz);
            mzBA(zi,pixi+(pixj-1)*m) = mzaxis(idz,pixj);
            % X((pixj-1)*np+(pixi-1)*n*np+zi) = z(idz);
            % Y((pixj-1)*np+(pixi-1)*n*np+zi) = mzaxis(idz,pixj);
        end
    end
end


zBlock = (reshape(zBA(round(np/2),:),m,n));
zBlock = zBlock(:,1);
xBlock = x;

X = zBA(:); % zBA
Y = mzBA(:); % mzBA

%% Gauss-Newton with LM
tol = 1e-3;
maxIte = 400;
muAlpha = 10;
muBeta = 10;
% muAlpha = 0; muBeta = 0;
beta0 = 1+(10.5)/2;
alpha0 = 0.1*freq^2; % dB/cm -> Np/m

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

estAClm = reshape(alphaArr/freq^2,[m,n]);
estBAlm = reshape(2*(betaArr-1),[m,n]);

fprintf('AC: %.3f +/- %.3f\n', mean(estAClm(:), 'omitnan'), ...
    std(estAClm(:), [] ,'omitnan'));
fprintf('B/A: %.2f +/- %.2f\n', mean(estBAlm(:), 'omitnan'), ...
    std(estBAlm(:), [] ,'omitnan'));

figure('Units','centimeters', 'Position',[5 5 10 5]), 
plot(loss, 'LineWidth',2)
xlim([5 length(loss)])
xlabel('Number of iterations')
ylabel('Loss')

figure; imagesc(x*1e3,z*1e3,estAClm); colorbar;
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


figure; imagesc(x*1e3,z*1e3,estBAlm); colorbar;
clim([5 10]);
title('B/A');
title(['B/A = ',num2str(median(estBAlm(:),'omitnan'),'%.1f')]);
font = 20;
axis image
colormap pink; colorbar;
set(gca,'fontsize',font)
xlabel('Lateral distance (mm)');
ylabel('Depth (mm)');
%title("\bf BoA map");
set(gca,'FontSize',font);
pause(0.5)


%% Getting mz
function [BAimage, BAparam] = BAestimator(media, param)

font = 10;
scale = param.factor;
noise_level = param.noise_level;
ACref = param.ACref;
ACsam = param.ACsam;
alpha_power_sam = param.APsam;
alpha_power_ref = param.APref;

load(media.samLP);
%rf = bf_planewave(sensor_data, fs, fnumber);
samLPfull = rf(:,:,1:end);
load(media.samHP);
%rf = bf_planewave(sensor_data, fs, fnumber);
samHPfull = rf(:,:,1:end); clear rf;
samLPfull = samLPfull + noise_level*randn(size(samLPfull));
samHPfull = samHPfull + noise_level*randn(size(samHPfull));

fs = param.fs;
%dt = 1/fs;
fL = param.f0 - 1.5; %1.5; %0.4*param.f0;
fH = param.f0 + 1.5; %8.5; %1.6*param.f0;
fnyquist = fs/2;

hfilter = fir1(param.order_filter,[(fL*1e6)/fnyquist (fH*1e6)/fnyquist]);

for ff = 1:size(samLPfull,3)
    for uu = 1:size(samHPfull,2)
        samLPfull(:,uu,ff) = conv(samLPfull(:,uu,ff),hfilter,'same');
        samHPfull(:,uu,ff) = conv(samHPfull(:,uu,ff),hfilter,'same');
    end
end


load(media.refLP);
%rf = bf_planewave(sensor_data, fs, fnumber);
refLPfull = rf(:,:,1); %refLPfull = scan_lines'; %refLP = refLP(:,col_last+1-width:col_last);
load(media.refHP);
%rf = bf_planewave(sensor_data, fs, fnumber);
refHPfull = rf(:,:,1); %refHPfull = scan_lines'; %refHP = refHP(:,col_last+1-width:col_last);
refLPfull = refLPfull + noise_level*randn(size(refLPfull));
refHPfull = refHPfull + noise_level*randn(size(refHPfull));


for uu = 1:size(refHPfull,2)
    refHPfull(:,uu) = conv(refHPfull(:,uu),hfilter,'same');
    refLPfull(:,uu) = conv(refLPfull(:,uu),hfilter,'same');
end


col_last = param.col_last;
overlap = param.overlap;
width = param.width;

n = param.n;
x_ini = x;
xBA = [];

while col_last <= 128

    samLP = samLPfull(:,col_last+1-width:col_last,:);
    samHP = samHPfull(:,col_last+1-width:col_last,:);

    refLP = refLPfull(:,col_last+1-width:col_last,:);
    refHP = refHPfull(:,col_last+1-width:col_last,:);

    refLP = refLPfull;
    refHP = refHPfull;

    %col_last - overlap
    xBA(end+1) = x(col_last - overlap);

    col_last = col_last + overlap;
    f0 = param.f0;
    dz = z(2)-z(1);
    L = param.size_wl * (c0/(f0*1e6)) / dz;


    if strcmp(param.instfreq,'no')
        f0_fund_sam = param.f0;
        f0_fund_ref = param.f0;
    end

    %% BoA estimation
    % x and z in cm
    %x = x*1e2;
    z = z*1e2;

    % z in meters again
    z = z(:)*1e-2;

    envLPsam = mean(abs(hilbert((samLP))),[2 3]);
    envHPsam = mean(abs(hilbert((samHP))),[2 3]);
    envLPref = mean(abs(hilbert((refLP))),[2 3]);
    envHPref = mean(abs(hilbert((refHP))),[2 3]);

    GAPsam = mean( (scale*envLPsam), 2) - mean( (envHPsam) ,2);
    GAPref = mean( (scale*envLPref), 2) - mean( (envHPref) ,2);

    alpha1 = ACsam;

    %f0_fund_sam = ifq_sam_ok;
    %f0_fund_sam = f0 - sigma_1st*sigma50pct^2 * f0^2 * ( alpha1.*(z*1e2) );
    ATT1MAPsam = 1./exp( mean( alpha1.* f0_fund_sam.^alpha_power_sam .*(z*1e2) , 2 ) );

    % f0_fund_ref = ifq_ref_ok;
    %f0_fund_ref = f0 - sigma_1st*sigma50pct^2 * f0^2 * ACref*(z*1e2);
    ATT1ref = 1./exp( ACref.*(f0_fund_ref.^alpha_power_ref) .*(z*1e2));

    P22sam = sqrt(abs(GAPsam).*envLPref);
    P22ref = sqrt(abs(GAPref).*envLPsam);



    P22sam_smooth = movmean(P22sam,L);
    P22ref_smooth = movmean(P22ref,L);


    alpha1_sam =  mean( ACsam .* f0_fund_sam.^alpha_power_sam , 2) *100;
    alpha1_ref =   ( ACref.* ( (f0_fund_ref).^alpha_power_ref)  )*100; % Np/m

    mz = ((1+0.5*param.BoAref))* (P22sam_smooth./P22ref_smooth) ...
        .* (1- ATT1ref.*ATT1ref)./(alpha1_ref.*(z));

    %ratio_beta2 = mz.* (alpha1_sam.*f0_fund_sam.^alpha_power_sam*(z))./(1- ATT1MAPsam.*ATT1MAPsam);

    ratio_beta = (P22sam_smooth./P22ref_smooth) ...
        .* (alpha1_sam.*z./(alpha1_ref.*z)) ...
        .* (1- ATT1ref.*ATT1ref)./(1- ATT1MAPsam.*ATT1MAPsam);

    BoAsam = ( ratio_beta*(1+0.5*param.BoAref) - 1 )*2;

    figure;
    plot(z, BoAsam,'r', 'LineWidth',3);
    xlim([5e-3 40e-3]); ylim([0 15]);
    %set(hfig,'units','normalized','outerposition',[0 0 1 1])
    %font = 20;
    set(gca,'fontsize',font)
    xlabel('\bfDepth (m)');
    %ylabel('\bfP_{22}^2');
    ylabel("\bf BoA ");
    set(gca,'FontSize',font);

    %end
    n = n+1;
    if col_last > 10
        disp(['n = ',num2str(n)]);
        %           keyboard
    end
    close all;


    BoAsam_matrix(:,n) = BoAsam;
    %keyboard
    mzmatrix(:,n) = mz;

end

%keyboard

BAimage = BoAsam_matrix;
BAparam.mz =  mzmatrix;
BAparam.z = z(:)*1e3;
%BAparam.x = x(:)*1e3;
BAparam.x = xBA(:)*1e3;

end


