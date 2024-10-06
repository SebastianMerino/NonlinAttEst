% My implementation of the depletion method with my own functions
% Used for testing several frequencies

clear; close all; clc;
addpath(genpath(pwd))

% baseDir = 'C:\Users\sebas\Documents\Data\Nonlinearity\rfdata';
% fileSam = "rf_fnum3_SCOMP5MHz_nc10_0p10f2_saminc400_doubleangle720";
% fileRef = "rf_fnum3_SCOMP5MHz_nc10_0p10f2_ref400_doubleangle720";
% baseDir = "C:\Users\sebas\Documents\Data\Nonlinearity\mySimulations";
% fileSam = "rf_baBack6_baInc9_att0p1";
% fileRef = "rf_ba8_att0p12_ref";
% samDir = 'C:\Users\sebas\Documents\Data\Nonlinearity\attInc\fn2';
% fileSam = "RFfn2_PWNE5MHz_samincBA6inc12_att0p1f2inc0p10_nc10_400kPa";
% refDir = "C:\Users\sebas\Documents\Data\Nonlinearity\homoAttMismatch\fn2";
% fileRef = "RFfn2_PWNE5MHz_refBA6_att0p10f2p0_400kPa";
samDir = 'C:\Users\sebas\Documents\Data\Nonlinearity\attInc\fn2';
fileSam = "RFfn2_PWNE5MHz_samincBA6inc12_att0p1f2inc0p10_nc10_400kPa";
refDir = 'C:\Users\sebas\Documents\Data\Nonlinearity\uniformBA_inclusionACS';
fileRef = "RFfn2_PWNE5MHz_samBA12_att0p1f2inc0p10_nc10_400kPa";

% Auxiliar variables
NptodB = 20*log10(exp(1));
radiusDisk = (9)*1e-3;
centerDepth = 22.5e-3;

% Hyperparameters
freqC = 5; 
v = 5; % scaling factor

% Known variables
betaR = 1 + 12/2;
alphaR = 0.1*freqC.^2/NptodB*100;
alphaS = 0.1*freqC.^2/NptodB*100;

%% Loading and cropping rf data
% Sample
sample = load(fullfile(samDir,fileSam));
z = sample.z';
x = sample.x;
fs = sample.fs;
rfL = sample.rf1(:,:,1);
rfH = sample.rf2(:,:,1);
% Reference
ref = load(fullfile(refDir,fileRef));
rfLR = ref.rf1(:,:,1);
rfHR = ref.rf2(:,:,1);


%% Filtering
freqTol = 0.5;
nCycles = 10; % Number of cycles of the initial filter

order = round(nCycles/freqC/1e6*fs);
PLfull = getFilteredPressure(rfL,fs,freqC*1e6,freqTol*1e6,order);
PHfull = getFilteredPressure(rfH,fs,freqC*1e6,freqTol*1e6,order);
PLRfull = getFilteredPressure(rfLR,fs,freqC*1e6,freqTol*1e6,order);
PHRfull = getFilteredPressure(rfHR,fs,freqC*1e6,freqTol*1e6,order);


% Plotting Bmode
Bmode = db(PLfull(z>0.5e-2,:));
Bmode = Bmode - max(Bmode(:));
figure, imagesc(x*100,z(z>0.5e-2)*100,Bmode, [-50 0])
hold on
rectangle('Position',[0-radiusDisk,centerDepth-radiusDisk,...
    2*radiusDisk,2*radiusDisk]*100, 'Curvature',1,...
    'EdgeColor','b', 'LineStyle','--', 'LineWidth',2)
hold off
axis image
colormap gray
xlabel('Lateral [cm]')
ylabel('Depth [cm]')

%% 
% Subsampling parameters
wl = 1540/freqC/1e6; % Mean central frequency
blockParams.blockSize = [10 10]*wl; 
blockParams.overlap = 0.8;
blockParams.zlim = [0.5; 5.5]/100;
blockParams.xlim = [-3; 3]/100;
[meanP,xP,zP] = getMeanBlock(cat(3,PLfull,PHfull,PLRfull,PHRfull),x,z,...
    blockParams);
PL = meanP(:,:,1);
PH = meanP(:,:,2);
PLR = meanP(:,:,3);
PHR = meanP(:,:,4);

% figure,imagesc(x*100,z*100,abs(v*PL-PH))
% axis image
% colorbar


%% Getting B/A
betaS = betaR*sqrt( abs(v*PL-PH)./abs(v*PLR-PHR) .*PLR./PL ).*...
    ( 1-exp(-2*alphaR*zP) )./( 1-exp(-2*alphaS*zP) ) *alphaS/alphaR;

% betaS = betaR*sqrt( abs(v*PL-PH)./abs(v*PLR-PHR) .*...
%     abs(v^3*PLR-PHR)./abs(v^3*PL-PH) ).*...
%     ( 1-exp(-2*alphaR*zP) )./( 1-exp(-2*alphaS*zP) ) *alphaS/alphaR;

BA = 2*(betaS-1);
% BA = 1./betaS;
cm = 100;
figure,
imagesc(xP*cm,zP*cm,BA, [5 14])
title(sprintf("B/A = %.2f \\pm %.2f",mean(BA(:)),std(BA(:))))
xlabel('Lateral [cm]')
ylabel('Depth [cm]')
axis image
colorbar
colormap pink
% hold on
% rectangle('Position',[0-radiusDisk,centerDepth-radiusDisk,...
%     2*radiusDisk,2*radiusDisk]*100, 'Curvature',1,...
%     'EdgeColor','b', 'LineStyle','--', 'LineWidth',2)
% hold off

%% Local B/A maps with regularization
[m,n]= size(BA);
muLocal = 0.01; %0.001;
tol = 1e-3;

dzP = zP(2)-zP(1);
izP = round(zP./dzP);

factorq = izP(1)./izP ;
estBAcum = BA - BA(1,:).*factorq; 
estBAcum = estBAcum(2:end,:);

P = sparse(tril(ones(m-1)));
P = P./izP(2:end);
P = kron(speye(n),P);
estBAinst = IRLS_TV(estBAcum(:),P,muLocal,m-1,n,tol,[],ones((m-1)*n,1));
estBAinst = reshape(estBAinst,m-1,n);

% estBAinst = 2*(1./estBAinst - 1);
baRange = [5,13];
figure; 
im = imagesc(xP*1e3,zP(2:end)*1e3,estBAinst); colorbar;
clim(baRange);
title('B/A');
axis image
colormap pink; colorbar;
xlabel('Lateral distance (mm)');
ylabel('Depth (mm)');
hold on
rectangle('Position',[0-radiusDisk,centerDepth-radiusDisk,...
2*radiusDisk,2*radiusDisk]*1000, 'Curvature',1,...
'EdgeColor','w', 'LineStyle','--', 'LineWidth',2)
hold off
% im.AlphaData = back;
pause(0.1)

[Xmesh,Zmesh] = meshgrid(xP,zP(2:end));
inc = Xmesh.^2 + (Zmesh-centerDepth).^2 < (radiusDisk-3e-3).^2;
back = Xmesh.^2 + (Zmesh-centerDepth).^2 > (radiusDisk+3e-3).^2;
fprintf('B/A inc: %.2f +/- %.2f\n', mean(estBAinst(inc),'omitnan'), ...
std(estBAinst(inc), [] ,'omitnan'));
fprintf('B/A back: %.2f +/- %.2f\n', mean(estBAinst(back),'omitnan'), ...
std(estBAinst(back), [] ,'omitnan'));
