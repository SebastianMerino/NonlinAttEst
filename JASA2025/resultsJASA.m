imPosition = [100 200 350 250];
lineWidth = 2;
fontSize = 12;
legendCell = {'GNTV','GNLM'};
legendLoc = 'northwest'; 
baseDir = 'Q:\smerino\Nonlinearity\resultsJASA';
baLim = [1 15];

%% Homogeneous data
tFile = fullfile(baseDir,'homoBA9','homoRef.xlsx');
T = readtable(tFile);
method = categorical(T.method);

baAdmm = T(method=='ADMM',:).baMean;
baStdAdmm = T(method=='ADMM',:).baStd;
baIus = T(method=='IUS',:).baMean;
baStdIus = T(method=='IUS',:).baStd;
alphaRef = T(method=='ADMM',:).alphaRef;
gammaRatio = (5.5./0.1)./(4./alphaRef);
xlimGamma = [gammaRatio(1)-0.3,gammaRatio(end)+0.3];

figure('Position',imPosition),
% errorbar(gammaRatio,baAdmm,baStdAdmm, 'LineWidth',lineWidth)
% hold on
% errorbar(gammaRatio,baIus,baStdIus, 'LineWidth',lineWidth)
errorbar(alphaRef,baAdmm,baStdAdmm, 'LineWidth',lineWidth)
hold on
errorbar(alphaRef,baIus,baStdIus, 'LineWidth',lineWidth)
hold off, grid on
yline(9,'k--')
title('Uniform media')
ylabel('B/A')
xlabel('\alpha_0 [dB/cm/MHz^2]')
legend(legendCell, 'Location',legendLoc)
% xlim(xlimGamma)
fontsize(fontSize,"points")
ylim([5 17])

acAdmm = T(method=='ADMM',:).acMean;
acStdAdmm = T(method=='ADMM',:).acStd;
acIus = T(method=='IUS',:).acMean;
acStdIus = T(method=='IUS',:).acStd;

figure('Position',imPosition),
% errorbar(gammaRatio,acAdmm,acStdAdmm, 'LineWidth',lineWidth)
% hold on
% errorbar(gammaRatio,acIus,acStdIus, 'LineWidth',lineWidth)
errorbar(alphaRef,acAdmm,acStdAdmm, 'LineWidth',lineWidth)
hold on
errorbar(alphaRef,acIus,acStdIus, 'LineWidth',lineWidth)
hold off, grid on
yline(0.1,'k--')
title('Uniform media')
ylabel('\alpha [dB/cm/MHz^2]')
xlabel('\alpha_0 [dB/cm/MHz^2]')
legend(legendCell, 'Location',legendLoc)
% xlim(xlimGamma)
fontsize(fontSize,"points")

%% Homogeneous data
tFile = fullfile(baseDir,'homoBA12','homoRef.xlsx');
T = readtable(tFile);
method = categorical(T.method);

baAdmm = T(method=='ADMM',:).baMean;
baStdAdmm = T(method=='ADMM',:).baStd;
baIus = T(method=='IUS',:).baMean;
baStdIus = T(method=='IUS',:).baStd;
alphaRef = T(method=='ADMM',:).alphaRef;
gammaRatio = (7./0.1)./(4./alphaRef);
xlimGamma = [gammaRatio(1)-0.3,gammaRatio(end)+0.3];

figure('Position',imPosition),
errorbar(gammaRatio,baAdmm,baStdAdmm, 'LineWidth',lineWidth)
hold on
errorbar(gammaRatio,baIus,baStdIus, 'LineWidth',lineWidth)
hold off, grid on
yline(12,'k--')
title('Uniform media')
ylabel('B/A')
xlabel('\Gamma_{s}/\Gamma_{r}')
% legend(legendCell, 'Location',legendLoc)
xlim(xlimGamma)
fontsize(fontSize,"points")
ylim([6 21])

acAdmm = T(method=='ADMM',:).acMean;
acStdAdmm = T(method=='ADMM',:).acStd;
acIus = T(method=='IUS',:).acMean;
acStdIus = T(method=='IUS',:).acStd;

figure('Position',imPosition),
errorbar(gammaRatio,acAdmm,acStdAdmm, 'LineWidth',lineWidth)
hold on
errorbar(gammaRatio,acIus,acStdIus, 'LineWidth',lineWidth)
hold off, grid on
yline(0.1,'k--')
title('Uniform media')
ylabel('\alpha [dB/cm/MHz^2]')
xlabel('\Gamma_{s}/\Gamma_{r}')
% legend(legendCell, 'Location',legendLoc)
xlim(xlimGamma)
fontsize(fontSize,"points")


%% BA6inc9, reference from uniformBA_inclusionACS
tFile = fullfile(baseDir,'ba6inc9Ref2p0','table.xlsx');
T = readtable(tFile);
method = categorical(T.method);
baIncAdmm = T(method=='ADMM',:).BaIncMean;
baIncStdAdmm = T(method=='ADMM',:).BaIncStd;
baIncIus = T(method=='IUS',:).BaIncMean;
baIncStdIus = T(method=='IUS',:).BaIncStd;
alphaInc = T(method=='ADMM',:).alphaInc;
gammaRatio = (5.5./alphaInc)/(4/0.1);
xlimGamma = [gammaRatio(end)-0.2,gammaRatio(1)+0.2];

figure('Position',imPosition),
errorbar(gammaRatio,baIncAdmm,baIncStdAdmm, 'LineWidth',lineWidth)
hold on
errorbar(gammaRatio,baIncIus,baIncStdIus, 'LineWidth',lineWidth)
hold off, grid on
yline(9,'k--')
title('Inclusion')
ylabel('B/A')
xlabel('\Gamma_{inc}/\Gamma_{bgnd}')
legend(legendCell, 'Location',legendLoc)
xlim(xlimGamma)
ylim(baLim)
fontsize(fontSize,"points")

acIncAdmm = T(method=='ADMM',:).AcIncMean;
acIncStdAdmm = T(method=='ADMM',:).AcIncStd;
acIncIus = T(method=='IUS',:).AcIncMean;
acIncStdIus = T(method=='IUS',:).AcIncStd;
alphaInc = T(method=='ADMM',:).alphaInc;

figure('Position',imPosition),
errorbar(gammaRatio,acIncAdmm-alphaInc,acIncStdAdmm, 'LineWidth',lineWidth)
hold on
errorbar(gammaRatio,acIncIus-alphaInc,acIncStdIus, 'LineWidth',lineWidth)
hold off, grid on
yline(0,'k--')
title('Inclusion')
ylabel('\Delta\alpha [dB/cm/MHz^2]')
xlabel('\Gamma_{inc}/\Gamma_{bgnd}')
legend(legendCell, 'Location',legendLoc)
xlim(xlimGamma)
ylim([-0.13 0.07])
fontsize(fontSize,"points")

baBackAdmm = T(method=='ADMM',:).BaBackMean;
baBackStdAdmm = T(method=='ADMM',:).BaBackStd;
baBackIus = T(method=='IUS',:).BaBackMean;
baBackStdIus = T(method=='IUS',:).BaBackStd;

figure('Position',imPosition),
errorbar(gammaRatio,baBackAdmm,baBackStdAdmm, 'LineWidth',lineWidth)
hold on
errorbar(gammaRatio,baBackIus,baBackStdIus, 'LineWidth',lineWidth)
hold off, grid on
yline(6,'k--')
title('Background')
ylabel('B/A')
xlabel('\Gamma_{inc}/\Gamma_{bgnd}')
% legend(legendCell, 'Location',legendLoc)
ylim(baLim)
xlim(xlimGamma)
fontsize(fontSize,"points")

acBackAdmm = T(method=='ADMM',:).AcBackMean;
acBackStdAdmm = T(method=='ADMM',:).AcBackStd;
acBackIus = T(method=='IUS',:).AcBackMean;
acBackStdIus = T(method=='IUS',:).AcBackStd;

figure('Position',imPosition),
errorbar(gammaRatio,acBackAdmm-0.1,acBackStdAdmm, 'LineWidth',lineWidth)
hold on
errorbar(gammaRatio,acBackIus-0.1,acBackStdIus, 'LineWidth',lineWidth)
hold off, grid on
yline(0,'k--')
title('Background')
ylabel('\Delta\alpha [dB/cm/MHz^2]')
xlabel('\Gamma_{inc}/\Gamma_{bgnd}')
% legend(legendCell, 'Location',legendLoc)
ylim([-0.13 0.07])
xlim(xlimGamma)
fontsize(fontSize,"points")

%% BA6inc12, reference from uniformBA_inclusionACS
tFile = fullfile(baseDir,'newSimulation','ba6inc12Ref2p0','table.xlsx');
T = readtable(tFile);
method = categorical(T.method);
baIncAdmm = T(method=='ADMM',:).BaIncMean;
baIncStdAdmm = T(method=='ADMM',:).BaIncStd;
baIncIus = T(method=='IUS',:).BaIncMean;
baIncStdIus = T(method=='IUS',:).BaIncStd;
alphaInc = T(method=='ADMM',:).alphaInc;
gammaRatio = (7./alphaInc)/(4/0.1);
xlimGamma = [gammaRatio(end)-0.2,gammaRatio(1)+0.2];

figure('Position',imPosition),
errorbar(gammaRatio,baIncAdmm,baIncStdAdmm, 'LineWidth',lineWidth)
hold on
errorbar(gammaRatio,baIncIus,baIncStdIus, 'LineWidth',lineWidth)
hold off, grid on
yline(12,'k--')
title('Inclusion')
ylabel('B/A')
xlabel('\Gamma_{inc}/\Gamma_{bgnd}')
% legend(legendCell, 'Location',legendLoc)
ylim(baLim)
xlim(xlimGamma)
fontsize(fontSize,"points")

acIncAdmm = T(method=='ADMM',:).AcIncMean;
acIncStdAdmm = T(method=='ADMM',:).AcIncStd;
acIncIus = T(method=='IUS',:).AcIncMean;
acIncStdIus = T(method=='IUS',:).AcIncStd;
alphaInc = T(method=='ADMM',:).alphaInc;

figure('Position',imPosition),
errorbar(gammaRatio,acIncAdmm-alphaInc,acIncStdAdmm, 'LineWidth',lineWidth)
hold on
errorbar(gammaRatio,acIncIus-alphaInc,acIncStdIus, 'LineWidth',lineWidth)
hold off, grid on
yline(0,'k--')
title('Inclusion')
ylabel('\Delta\alpha [dB/cm/MHz^2]')
xlabel('\Gamma_{inc}/\Gamma_{bgnd}')
% legend(legendCell, 'Location',legendLoc)
ylim([-0.13 0.07])
xlim(xlimGamma)
fontsize(fontSize,"points")

baBackAdmm = T(method=='ADMM',:).BaBackMean;
baBackStdAdmm = T(method=='ADMM',:).BaBackStd;
baBackIus = T(method=='IUS',:).BaBackMean;
baBackStdIus = T(method=='IUS',:).BaBackStd;

figure('Position',imPosition),
errorbar(gammaRatio,baBackAdmm,baBackStdAdmm, 'LineWidth',lineWidth)
hold on
errorbar(gammaRatio,baBackIus,baBackStdIus, 'LineWidth',lineWidth)
hold off, grid on
yline(6,'k--')
title('Background')
ylabel('B/A')
xlabel('\Gamma_{inc}/\Gamma_{bgnd}')
% legend(legendCell, 'Location',legendLoc)
ylim(baLim)
xlim(xlimGamma)
fontsize(fontSize,"points")

acBackAdmm = T(method=='ADMM',:).AcBackMean;
acBackStdAdmm = T(method=='ADMM',:).AcBackStd;
acBackIus = T(method=='IUS',:).AcBackMean;
acBackStdIus = T(method=='IUS',:).AcBackStd;

figure('Position',imPosition),
errorbar(gammaRatio,acBackAdmm-0.1,acBackStdAdmm, 'LineWidth',lineWidth)
hold on
errorbar(gammaRatio,acBackIus-0.1,acBackStdIus, 'LineWidth',lineWidth)
hold off, grid on
yline(0,'k--')
title('Background')
ylabel('\Delta\alpha [dB/cm/MHz^2]')
xlabel('\Gamma_{inc}/\Gamma_{bgnd}')
% legend(legendCell, 'Location',legendLoc)
ylim([-0.13 0.07])
xlim(xlimGamma)
fontsize(fontSize,"points")

%%
outDir = 'Q:\smerino\Nonlinearity\resultsJASA\plotsRef1p2';
mkdir(outDir)
save_all_figures_to_directory(outDir,'fig','svg')
close all