% My implementation of the depletion method with my own functions
% Used for testing several frequencies

clear; close all; clc;
addpath(genpath(pwd))

freqC = 5;
samDir = 'C:\Users\sebas\Documents\Data\Nonlinearity\attInc\fn2';
resultsDir = fullfile(samDir,"results","ba6inc9DepletionRefBack");
[~,~,~] = mkdir(resultsDir);

refDir = "C:\Users\sebas\Documents\Data\Nonlinearity\homoAttMismatch\fn2";
fileRef = "RFfn2_PWNE5MHz_refBA6_att0p10f2p0_400kPa";
% refDir = 'C:\Users\sebas\Documents\Data\Nonlinearity\uniformBA_inclusionACS';
% fileRef = "RFfn2_PWNE"+freqC+"MHz_samBA12_att0p1f2inc0p10_nc10_400kPa";

% Auxiliar variables
NptodB = 20*log10(exp(1));
radiusDisk = (9)*1e-3;
centerDepth = 22.5e-3;

% Hyperparameters
v = 5; % scaling factor

% Known variables
betaR = 1 + 6/2;
alphaR = 0.1*freqC.^2/NptodB*100;
alphaS = 0.1*freqC.^2/NptodB*100;
attRange = [0.08,0.20];
%%
iSim = 2;
alphaIncVec = 0.08:0.02:0.22;

%%
for iSim = 1:length(alphaIncVec)
alphaInc = alphaIncVec(iSim);
alphaStr = num2str(round(alphaInc*100),"%02d");

fileSam = "RFfn2_PWNE"+freqC+"MHz_samincBA6inc9_att0p1f2inc0p"+alphaStr+ ...
    "_nc10_400kPa";


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
Bmode = db(PLfull(z>0.5e-2 & z<5.5e-2,:));
Bmode = Bmode - max(Bmode(:));
figure, imagesc(x*100,z(z>0.5e-2 & z<5.5e-2)*100,Bmode, [-50 0])
hold on
rectangle('Position',[0-radiusDisk,centerDepth-radiusDisk,...
    2*radiusDisk,2*radiusDisk]*100, 'Curvature',1,...
    'EdgeColor','b', 'LineStyle','--', 'LineWidth',2)
hold off
axis image
colormap gray
xlabel('Lateral [cm]')
ylabel('Depth [cm]')
title('Bmode')
%%
% Subsampling parameters
wl = 1540/freqC/1e6; % Mean central frequency
blockParams.blockSize = [15 15]*wl;
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

%% Alpha maps
alphaIncNpm = alphaInc*freqC.^2/NptodB*100;
[X,Z] = meshgrid(xP,zP);

inc = X.^2 + (Z-centerDepth).^2 <= (radiusDisk).^2;
[m,n] = size(PL);
alphaL = ones(size(PL))*alphaS;
alphaL(inc) = alphaIncNpm;
izBlock = 1:length(zP);
P = sparse(tril(ones(m)));
P = P./izBlock';
alphaC =  P*(alphaL);

figure,
% tiledlayout(1,2)
% nexttile,
imagesc(xP*1e3,zP*1e3,alphaC/freqC.^2*NptodB/100, attRange)
title('Alpha cumulative')
axis image
colormap turbo; colorbar;
xlabel('Lateral distance (mm)');
ylabel('Depth (mm)');
% nexttile,
% imagesc(xP,zP,alphaSZ)
% axis image

alphaSZ = alphaC.*Z;
alphaRZ = alphaR.*Z;

%% Getting B/A
% betaS = betaR*sqrt( abs(v*PL-PH)./abs(v*PLR-PHR) .*PLR./PL ).*...
%     ( 1-exp(-2*alphaR*zP) )./( 1-exp(-2*alphaS*zP) ) *alphaS/alphaR;


betaS = betaR*sqrt( abs(v*PL-PH)./abs(v*PLR-PHR) .*PLR./PL ).*...
    ( 1-exp(-2*alphaRZ) )./( 1-exp(-2*alphaSZ) ) .*alphaSZ./alphaRZ;

BA = 2*(betaS-1);
cm = 100;
% figure,
% imagesc(xP*cm,zP*cm,BA, [5 14])
% title(sprintf("B/A = %.2f \\pm %.2f",mean(BA(:)),std(BA(:))))
% xlabel('Lateral [cm]')
% ylabel('Depth [cm]')
% axis image
% colorbar
% colormap pink
% hold on
% rectangle('Position',[0-radiusDisk,centerDepth-radiusDisk,...
%     2*radiusDisk,2*radiusDisk]*100, 'Curvature',1,...
%     'EdgeColor','b', 'LineStyle','--', 'LineWidth',2)
% hold off

%% Local B/A maps with regularization
[m,n]= size(BA);
muLocal = 0.1; %0.001;
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

baRange = [5,13];

%% ROI and metrics
%  Masks
xBm = x; zBm = z;
[Xq, Zq] = meshgrid(xBm,zBm);
L = radiusDisk; cx = radiusDisk*1.5;
inc = maskRect(xBm, zBm, 0, centerDepth, L, L);
back = maskRect(xBm, zBm, -cx, centerDepth, L/2, L) | ...
    maskRect(xBm, zBm, cx, centerDepth, L/2, L);

% Interp
[Xmesh,Zmesh] = meshgrid(xP,zP(2:end));
baInterp = interp2(Xmesh,Zmesh,estBAinst,Xq,Zq);

% Metrics
metrics(iSim) = getMetrics(baInterp,inc,back,'Dep',alphaInc);

%% Plots
figure; imagesc(xP*1e3,zP(2:end)*1e3,estBAinst); colorbar;
clim(baRange);
title('Depletion, B/A');
axis image
colormap pink; colorbar;
xlabel('Lateral distance (mm)');
ylabel('Depth (mm)');
hold on
rectangle('Position',[0-radiusDisk,centerDepth-radiusDisk,...
    2*radiusDisk,2*radiusDisk]*1000, 'Curvature',1,...
    'EdgeColor','w', 'LineStyle','--', 'LineWidth',2)
rectangle('Position',[-L/2,centerDepth-L/2,...
    L,L]*1000, 'EdgeColor','k', 'LineStyle','--', 'LineWidth',2)
rectangle('Position',[cx-L/4,centerDepth-L/2,...
    L/2,L]*1000, 'EdgeColor','k', 'LineStyle','--', 'LineWidth',2)
hold off
rectangle('Position',[-cx-L/4,centerDepth-L/2,...
    L/2,L]*1000, 'EdgeColor','k', 'LineStyle','--', 'LineWidth',2)
hold off
pause(0.1)
%%
    save_all_figures_to_directory(resultsDir,...
        char("sim"+iSim+"fig"))
    close all

end

T = [struct2table(metrics)];
writetable(T,fullfile(resultsDir,'ba6inc9Depletion.xlsx'))

%%
function metrics = getMetrics(BA,inc,back,method,alphaInc)
metrics.BaIncMean = mean(BA(inc), 'omitnan');
metrics.BaIncStd = std(BA(inc), [], 'omitnan');
metrics.BaBackMean = mean(BA(back), 'omitnan');
metrics.BaBackStd = std(BA(back), [], 'omitnan');
metrics.method = method;
metrics.alphaInc = alphaInc;
end
