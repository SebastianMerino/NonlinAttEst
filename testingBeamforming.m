baseDir = 'C:\Users\smerino.C084288\Documents\MATLAB\BA_AC_joint\rfdata';
fileSam = 'rf_prebf_BaAttRef.mat';
load(fullfile(baseDir,fileSam))

focalNumber = 3;
dz = z(2) - z(1);
pitch = x(2) - x(1);


rf1 = BFangle(rf1,length(x),fs,1540, pitch,'bh',focalNumber,0);
rf2 = BFangle(rf2,length(x),fs,1540, pitch,'bh',focalNumber,0);


%%
save(fullfile(baseDir,'rf_ref_ba8_att0p12.mat'),...
    'rf1','rf2','x','z','fs');
%%
offset = 0.15; % deadband in cm
offset = round(offset/100/dz)+1;
rf1 = rf1(offset:end,:);
z = z(offset:end);
bmode1 = db(hilbert(rf1));
bmode1 = bmode1 - max(bmode1(:));

figure,
imagesc(x*100,z*100,bmode1, [-40 0])
colormap gray
axis image
