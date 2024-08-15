BaLimits = [7 11];

resultDir = fullfile('results','muBdivmuA1');

ini1 = load(fullfile(resultDir,'samBA9_attIni0p08'));
ini2 = load(fullfile(resultDir,'samBA9_attIni0p10'));
ini3 = load(fullfile(resultDir,'samBA9_attIni0p12'));

figure
errorbar(ini1.mu_vector,ini1.mBA,ini1.sBA, 'Color',[1 0.5 0], 'LineWidth',3)
set(gca, 'XScale', 'log')
hold on;
errorbar(ini2.mu_vector,ini2.mBA,ini2.sBA, 'Color',[1 0 0], 'LineWidth',3)
errorbar(ini3.mu_vector,ini3.mBA,ini3.sBA, 'Color',[1 0 0.5], 'LineWidth',3)
plot(mu_vector,9*ones(size(mu_vector)),'k--','LineWidth',3)
xlabel('\mu ')
title('B/A')
ylim(BaLimits), 
xlim([1e-4 1e5])
legend('\alpha_{ini} = 0.08','\alpha_{ini} = 0.10','\alpha_{ini} = 0.12',...
    'Location','northoutside')

%%
resultDir = fullfile('results','muBdivmuA10');

ini1 = load(fullfile(resultDir,'samBA9_attIni0p08'));
ini2 = load(fullfile(resultDir,'samBA9_attIni0p10'));
ini3 = load(fullfile(resultDir,'samBA9_attIni0p12'));

figure
errorbar(ini1.mu_vector,ini1.mBA,ini1.sBA, 'Color',[1 0.5 0], 'LineWidth',3)
set(gca, 'XScale', 'log')
hold on;
errorbar(ini2.mu_vector,ini2.mBA,ini2.sBA, 'Color',[1 0 0], 'LineWidth',3)
errorbar(ini3.mu_vector,ini3.mBA,ini3.sBA, 'Color',[1 0 0.5], 'LineWidth',3)
plot(mu_vector,9*ones(size(mu_vector)),'k--','LineWidth',3)
xlabel('\mu ')
title('B/A')
ylim(BaLimits), 
xlim([1e-4 1e5])
legend('\alpha_{ini} = 0.08','\alpha_{ini} = 0.10','\alpha_{ini} = 0.12',...
    'Location','northoutside')

%%
resultDir = fullfile('results','muBdivmuA100');

ini1 = load(fullfile(resultDir,'samBA9_attIni0p08'));
ini2 = load(fullfile(resultDir,'samBA9_attIni0p10'));
ini3 = load(fullfile(resultDir,'samBA9_attIni0p12'));

figure
errorbar(ini1.mu_vector,ini1.mBA,ini1.sBA, 'Color',[1 0.5 0], 'LineWidth',3)
set(gca, 'XScale', 'log')
hold on;
errorbar(ini2.mu_vector,ini2.mBA,ini2.sBA, 'Color',[1 0 0], 'LineWidth',3)
errorbar(ini3.mu_vector,ini3.mBA,ini3.sBA, 'Color',[1 0 0.5], 'LineWidth',3)
plot(mu_vector,9*ones(size(mu_vector)),'k--','LineWidth',3)
xlabel('\mu ')
title('B/A')
ylim(BaLimits),
xlim([1e-4 1e5])
legend('\alpha_{ini} = 0.08','\alpha_{ini} = 0.10','\alpha_{ini} = 0.12',...
    'Location','northoutside')
