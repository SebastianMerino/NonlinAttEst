filename = 'rf_prebf_BaAttInc1.mat';
load(filename)

focalNumber = 1;
dz = z(2) - z(1);
offset = 0.15; % deadband in cm
offset = round(offset/100/dz)+1;

rf1 = BFangle(rf1,length(x),fs,1540, pitch,'bh',focalNumber,0);
rf2 = BFangle(rf2,length(x),fs,1540, pitch,'bh',focalNumber,0);

rf1 = rf1(offset:end,:);
z = z(offset:end);
bmode1 = db(hilbert(rf1));
bmode1 = bmode1 - max(bmode1(:));

figure,
imagesc(x*100,z*100,bmode1, [-40 0])
colormap gray
axis image

%%
save(['rf_',filename(9:end)],...
    'rf1','rf2','x','z','fs');
