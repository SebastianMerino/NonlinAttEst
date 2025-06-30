clear, clc
dataDir = 'Q:\smerino\Nonlinearity\AC_UiX_new';
resultsDir = fullfile(dataDir,"bf");
mkdir(resultsDir)
filesLP = dir(fullfile(dataDir,"PWFOC*5MHz*80kPa.mat"));
fnumber = 2;
freq0 = 5e6;
nCycles = 10;

%%
for iFile = 1:length(filesLP)
    clear rf1 rf2
    file_LP = filesLP(iFile).name;
    file_HP = [filesLP(iFile).name(1:end-9),'400kPa.mat'];
    file_out = ['RFfn',num2str(fnumber),'_',file_HP];
    disp(file_out)

    tic
    load(fullfile(dataDir,file_LP))
    t0 = nCycles/freq0/2;
    z = (1:length(rf_prebf))*(1/fs)*c0/2;
    x = (1:size(rf_prebf,2))*pitch; x = x - mean(x);
    parfor rr = 1:size(rf_prebf,3)
        rf1(:,:,rr) = bfPlaneWaveAcq(rf_prebf(:,:,rr),x,z,t0,c0,fs,fnumber);
    end
    toc

    tic
    load(fullfile(dataDir,file_HP))
    parfor rr = 1:size(rf_prebf,3)
        rf2(:,:,rr) = bfPlaneWaveAcq(rf_prebf(:,:,rr),x,z,t0,c0,fs,fnumber);
    end
    toc
    
    zValid = z>1e-3 & z< 5.4e-2; 
    bMode = getBmode(rf1(zValid,:,1),fs,[1e6,10e6]);
    figure,
    imagesc(x*100,z(zValid)*100,bMode, [-60 0])
    colormap gray
    axis image
    pause(0.1)

    save(fullfile(resultsDir,file_out),'c0','fs','rf1','rf2','x','z');

end

clear rf1 rf2


%% Check
% for aa=1:4
%     zValid = z>1e-3 & z< 5.4e-2; 
%     bMode(:,:,aa) = getBmode(rf1(zValid,:,aa),fs,[1e6,10e6]);
%     figure,
%     imagesc(x*100,z(zValid)*100,bMode(:,:,aa), [-60 0])
%     colormap gray
%     axis image
%     pause(0.1)
% end
% cuec = medfilt2( mean(bMode,3),[420,7], 'symmetric');
% figure,
% imagesc(x*100,z(zValid)*100,cuec)
% colormap parula
% axis image
% pause(0.1)
