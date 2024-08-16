baseDir = 'C:\Users\smerino.C084288\Documents\MATLAB\BA_AC_joint\rfdata';
fileSam = 'rf_prebf_BaAttRef.mat';
load(fullfile(baseDir,fileSam))

focalNumber = 2;
dz = z(2) - z(1);
pitch = x(2) - x(1);
offset = 0.15; % deadband in cm
offset = round(offset/100/dz)+1;

rf1 = BFangle(rf1,length(x),fs,1540, pitch,'bh',focalNumber,0);
rf2 = BFangle(rf2,length(x),fs,1540, pitch,'bh',focalNumber,0);

save(fullfile(baseDir,'rf_ba8_att0p12_ref.mat'),...
    'rf1','rf2','x','z','fs');

rf1 = rf1(offset:end,:);
z = z(offset:end);
bmode1 = db(hilbert(rf1));
bmode1 = bmode1 - max(bmode1(:));

figure,
imagesc(x*100,z*100,bmode1, [-40 0])
colormap gray
axis image
%%
baseDir = 'C:\Users\smerino.C084288\Documents\MATLAB\BA_AC_joint\rfdata';
fileSam = 'rf_prebf_BaAttInc1.mat';
load(fullfile(baseDir,fileSam))

rf1 = BFangle(rf1,length(x),fs,1540, pitch,'bh',focalNumber,0);
rf2 = BFangle(rf2,length(x),fs,1540, pitch,'bh',focalNumber,0);

save(fullfile(baseDir,'rf_baBack6_baInc9_att0p1.mat'),...
    'rf1','rf2','x','z','fs');

rf1 = rf1(offset:end,:);
z = z(offset:end);
bmode1 = db(hilbert(rf1));
bmode1 = bmode1 - max(bmode1(:));

figure,
imagesc(x*100,z*100,bmode1, [-40 0])
colormap gray
axis image
%%
baseDir = 'C:\Users\smerino.C084288\Documents\MATLAB\BA_AC_joint\rfdata';
fileSam = 'rf_prebf_BaAttInc2.mat';
load(fullfile(baseDir,fileSam))

rf1 = BFangle(rf1,length(x),fs,1540, pitch,'bh',focalNumber,0);
rf2 = BFangle(rf2,length(x),fs,1540, pitch,'bh',focalNumber,0);

save(fullfile(baseDir,'rf_ba9_attBack0p1_attInc0p18.mat'),...
    'rf1','rf2','x','z','fs');

rf1 = rf1(offset:end,:);
z = z(offset:end);
bmode1 = db(hilbert(rf1));
bmode1 = bmode1 - max(bmode1(:));

figure,
imagesc(x*100,z*100,bmode1, [-40 0])
colormap gray
axis image

%%
baseDir = 'C:\Users\smerino.C084288\Documents\MATLAB\BA_AC_joint\rfdata';
fileSam = 'rf_prebf_BaAttInc3.mat';
load(fullfile(baseDir,fileSam))

rf1 = BFangle(rf1,length(x),fs,1540, pitch,'bh',focalNumber,0);
rf2 = BFangle(rf2,length(x),fs,1540, pitch,'bh',focalNumber,0);

save(fullfile(baseDir,'rf_baBack6_baInc9_attBack0p1_attInc0p18.mat'),...
    'rf1','rf2','x','z','fs');

rf1 = rf1(offset:end,:);
z = z(offset:end);
bmode1 = db(hilbert(rf1));
bmode1 = bmode1 - max(bmode1(:));

figure,
imagesc(x*100,z*100,bmode1, [-40 0])
colormap gray
axis image