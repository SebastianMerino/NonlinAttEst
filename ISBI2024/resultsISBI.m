imPosition = [100 200 250 250];
lineWidth = 2;
fontSize = 12;
legendCell = {'GNTV','GNLM'};
legendLoc = 'northwest'; 
baseDir = 'Q:\smerino\Nonlinearity\resultsISBI';
baLim = [1 15];

%% BA6inc12, reference from uniformBA_inclusionACS
tFile = fullfile(baseDir,'ba6Inc12','table.xlsx');
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
legend(legendCell, 'Location',legendLoc)
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
legend(legendCell, 'Location',legendLoc)
ylim([-0.13 0.07])
xlim(xlimGamma)
fontsize(fontSize,"points")

% %% BA6inc12, reference from uniformBA_inclusionACS
% tFile = fullfile(baseDir,'ba6Inc12','table.xlsx');
% T = readtable(tFile);
% method = categorical(T.method);
% baIncAdmm = T(method=='ADMM',:).BaIncMean;
% baIncStdAdmm = T(method=='ADMM',:).BaIncStd;
% baIncIus = T(method=='IUS',:).BaIncMean;
% baIncStdIus = T(method=='IUS',:).BaIncStd;
% alphaInc = T(method=='ADMM',:).alphaInc;
% xlimGamma = [alphaInc(1)-0.01,alphaInc(end)+0.01];
% 
% figure('Position',imPosition),
% errorbar(alphaInc,baIncAdmm,baIncStdAdmm, 'LineWidth',lineWidth)
% hold on
% errorbar(alphaInc,baIncIus,baIncStdIus, 'LineWidth',lineWidth)
% hold off, grid on
% yline(12,'k--')
% title('Inclusion')
% ylabel('B/A')
% xlabel('\alpha_{inc} [dB/cm/MHz^2]')
% % legend(legendCell, 'Location',legendLoc)
% ylim(baLim)
% xlim(xlimGamma)
% fontsize(fontSize,"points")
% 
% acIncAdmm = T(method=='ADMM',:).AcIncMean;
% acIncStdAdmm = T(method=='ADMM',:).AcIncStd;
% acIncIus = T(method=='IUS',:).AcIncMean;
% acIncStdIus = T(method=='IUS',:).AcIncStd;
% alphaInc = T(method=='ADMM',:).alphaInc;
% 
% figure('Position',imPosition),
% errorbar(alphaInc,acIncAdmm,acIncStdAdmm, 'LineWidth',lineWidth)
% hold on
% errorbar(alphaInc,acIncIus,acIncStdIus, 'LineWidth',lineWidth)
% plot(alphaInc,alphaInc,'k--')
% hold off, grid on
% title('Inclusion')
% ylabel('\alpha [dB/cm/MHz^2]')
% xlabel('\alpha_{inc} [dB/cm/MHz^2]')
% % legend(legendCell, 'Location',legendLoc)
% ylim([0.05 0.25])
% xlim(xlimGamma)
% fontsize(fontSize,"points")
% 
% baBackAdmm = T(method=='ADMM',:).BaBackMean;
% baBackStdAdmm = T(method=='ADMM',:).BaBackStd;
% baBackIus = T(method=='IUS',:).BaBackMean;
% baBackStdIus = T(method=='IUS',:).BaBackStd;
% 
% figure('Position',imPosition),
% errorbar(alphaInc,baBackAdmm,baBackStdAdmm, 'LineWidth',lineWidth)
% hold on
% errorbar(alphaInc,baBackIus,baBackStdIus, 'LineWidth',lineWidth)
% hold off, grid on
% yline(6,'k--')
% title('Background')
% ylabel('B/A')
% xlabel('\alpha_{inc} [dB/cm/MHz^2]')
% ylim(baLim)
% xlim(xlimGamma)
% fontsize(fontSize,"points")
% 
% acBackAdmm = T(method=='ADMM',:).AcBackMean;
% acBackStdAdmm = T(method=='ADMM',:).AcBackStd;
% acBackIus = T(method=='IUS',:).AcBackMean;
% acBackStdIus = T(method=='IUS',:).AcBackStd;
% 
% figure('Position',imPosition),
% errorbar(alphaInc,acBackAdmm,acBackStdAdmm, 'LineWidth',lineWidth)
% hold on
% errorbar(alphaInc,acBackIus,acBackStdIus, 'LineWidth',lineWidth)
% hold off, grid on
% yline(0.1,'k--')
% title('Background')
% ylabel('\alpha [dB/cm/MHz^2]')
% xlabel('\alpha_{inc} [dB/cm/MHz^2]')
% % legend(legendCell, 'Location',legendLoc)
% ylim([0.05 0.25])
% xlim(xlimGamma)
% fontsize(fontSize,"points")


%%
outDir = 'Q:\smerino\Nonlinearity\resultsISBI\plots';
mkdir(outDir)
save_all_figures_to_directory(outDir,'fig','svg')
close all