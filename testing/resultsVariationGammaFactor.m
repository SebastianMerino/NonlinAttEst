imPosition = [100 200 350 250];
lineWidth = 2;
fontSize = 12;
legendCell = {'GNTV'};
legendLoc = 'northwest'; 
baseDir = 'P:\smerino\NonlinAttESt\data\ref_att1p2';
baLim = [6 14];

%% Homogeneous data
tFile = fullfile(baseDir,'results_BA10_att0p60f1p2_gamma1p2_v4p5','homoRef.xlsx');
T = readtable(tFile);
method = categorical(T.method);

baAdmm = T(method=='ADMM',:).baMean;
baStdAdmm = T(method=='ADMM',:).baStd;
baIus = T(method=='IUS',:).baMean;
baStdIus = T(method=='IUS',:).baStd;
alphaRef = T(method=='ADMM',:).alphaRef;
gammaRatio = (6./(0.6*5^1.2))./(4./(alphaRef*5^2));
% xlimGamma = [gammaRatio(1)-0.3,gammaRatio(end)+0.3];

figure('Position',imPosition),
errorbar(gammaRatio,baAdmm,baStdAdmm, 'LineWidth',lineWidth)
hold on
% errorbar(gammaRatio,baIus,baStdIus, 'LineWidth',lineWidth)
hold off, grid on
yline(10,'k--')
title('Uniform media')
ylabel('B/A')
% xlabel('\alpha_{ref}')
xlabel('\Gamma_{sam}/\Gamma_{ref}')
legend(legendCell, 'Location',legendLoc)
% xlim(xlimGamma)
fontsize(fontSize,"points")
ylim(baLim)