clear, clc
dataDir = 'Q:\smerino\Nonlinearity\attIncNonQuadratic';
filesLP = dir(fullfile(dataDir,"*_80kPa.mat"));
fnumber = 2;
element_pitch = 0.3e-3;
%iFile = 1;
%%
for iFile = 1:length(filesLP)
    file_LP = filesLP(iFile).name;
    file_HP = [filesLP(iFile).name(1:end-9),'400kPa.mat'];
    file_out = ['RFfn',num2str(fnumber),'_',file_HP];
    disp(file_out)

    tic
    load(fullfile(dataDir,file_LP))
    z = (1:length(rf_prebf))*(1/fs)*c0/2;
    x = (1:size(rf_prebf,2))*element_pitch; x = x - mean(x);
    parfor rr = 1:size(rf_prebf,3)
        rf1(:,:,rr) = bf_planewave((rf_prebf(:,:,rr))',fs,fnumber,c0);
    end
    toc

    tic
    load(fullfile(dataDir,file_HP))
    parfor rr = 1:size(rf_prebf,3)
        rf2(:,:,rr) = bf_planewave((rf_prebf(:,:,rr))',fs,fnumber,c0);
    end
    toc

    save(fullfile(dataDir,file_out),'c0','fs','rf1','rf2','x','z');

    bMode = db(hilbert(rf1(:,:,1)));
    bMode = bMode - max(bMode(:));
    figure,
    imagesc(x,z,bMode, [-90 -30])
    colormap gray
    axis image
end

clear rf1 rf2


%%
offset = 500;
bmode1 = db(hilbert(rf1(offset:end,:,1)));
bmode1 = bmode1 - max(bmode1(:));
bmode2 = db(hilbert(rf2(offset:end,:,1)));
bmode2 = bmode2 - max(bmode2(:));

figure,
imagesc(x,z(offset:end),bmode1, [-50 0])
axis image
colorbar
colormap gray

figure,
imagesc(x,z(offset:end),bmode2, [-50 0])
axis image
colorbar
colormap gray