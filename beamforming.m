setup,

baseDir = 'C:\Users\smerino.C084288\Documents\Datasets\Nonlinearity';
fileSam = 'rf_prebf_BaAttIncFreq5.mat';
load(fullfile(baseDir,fileSam))

focalNumber = 2;
dz = z(2) - z(1);
pitch = x(2) - x(1);
offset = 0.15; % deadband in cm
offset = round(offset/100/dz)+1;

rf1 = BFangle(rf1,length(x),fs,1540, pitch,'bh',focalNumber,0);
rf2 = BFangle(rf2,length(x),fs,1540, pitch,'bh',focalNumber,0);

save(fullfile(baseDir,'rf_BaAttIncFreq5.mat'),...
    'rf1','rf2','x','z','fs');

rf1 = rf1(offset:end,:);
z = z(offset:end);
bmode1 = db(hilbert(rf1));
bmode1 = bmode1 - max(bmode1(:));

figure,
imagesc(x*100,z*100,bmode1, [-50 0])
colormap gray
rectangle('Position',[-0.8,2.5-0.8,1.6,1.6], 'Curvature',1, ...
    'EdgeColor','w', 'LineStyle','--', 'LineWidth',1)
axis image
%%
fileSam = 'rf_prebf_BaAttIncFreq6.mat';
load(fullfile(baseDir,fileSam))

rf1 = BFangle(rf1,length(x),fs,1540, pitch,'bh',focalNumber,0);
rf2 = BFangle(rf2,length(x),fs,1540, pitch,'bh',focalNumber,0);

save(fullfile(baseDir,'rf_BaAttIncFreq6.mat'),...
    'rf1','rf2','x','z','fs');

rf1 = rf1(offset:end,:);
z = z(offset:end);
bmode1 = db(hilbert(rf1));
bmode1 = bmode1 - max(bmode1(:));

figure,
imagesc(x*100,z*100,bmode1, [-50 0])
colormap gray
rectangle('Position',[-0.8,2.5-0.8,1.6,1.6], 'Curvature',1, ...
    'EdgeColor','w', 'LineStyle','--', 'LineWidth',1)
axis image
%%
fileSam = 'rf_prebf_BaAttIncFreq7.mat';
load(fullfile(baseDir,fileSam))

rf1 = BFangle(rf1,length(x),fs,1540, pitch,'bh',focalNumber,0);
rf2 = BFangle(rf2,length(x),fs,1540, pitch,'bh',focalNumber,0);

save(fullfile(baseDir,'rf_BaAttIncFreq7.mat'),...
    'rf1','rf2','x','z','fs');

rf1 = rf1(offset:end,:);
z = z(offset:end);
bmode1 = db(hilbert(rf1));
bmode1 = bmode1 - max(bmode1(:));

figure,
imagesc(x*100,z*100,bmode1, [-50 0])
colormap gray
rectangle('Position',[-0.8,2.5-0.8,1.6,1.6], 'Curvature',1, ...
    'EdgeColor','w', 'LineStyle','--', 'LineWidth',1)
axis image
%%
% fileSam = 'rf_prebf_BaAttInc3.mat';
% load(fullfile(baseDir,fileSam))
% 
% rf1 = BFangle(rf1,length(x),fs,1540, pitch,'bh',focalNumber,0);
% rf2 = BFangle(rf2,length(x),fs,1540, pitch,'bh',focalNumber,0);
% 
% save(fullfile(baseDir,'rf_baBack6_baInc9_attBack0p1_attInc0p18.mat'),...
%     'rf1','rf2','x','z','fs');
% 
% rf1 = rf1(offset:end,:);
% z = z(offset:end);
% bmode1 = db(hilbert(rf1));
% bmode1 = bmode1 - max(bmode1(:));
% 
% figure,
% imagesc(x*100,z*100,bmode1, [-40 0])
% colormap gray
% axis image