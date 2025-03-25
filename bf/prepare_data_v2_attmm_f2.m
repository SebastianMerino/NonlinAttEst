clear; close all; clc;
dataDir = 'D:\NonlinAttESt\data\nonUniformSimulation\ref';

x = linspace(0,0.0384,128); x = x-mean(x);

freqString = {'4p5','5p5'};
for ba = 6 % 9:12
    for alpha = 10 % :2:22
        for iFreq = 1
            freq = freqString{iFreq};
            % file_LP = ['PWNE',freq,'MHz_samincBA6inc',num2str(ba),...
            %     '_att0p1f2inc0p',num2str(alpha,'%02d'),'_nc10_80kPa'];
            % file_HP = ['PWNE',freq,'MHz_samincBA6inc',num2str(ba),...
            %     '_att0p1f2inc0p',num2str(alpha,'%02d'),'_nc10_400kPa'];
            file_LP = ['PWNE',freq,'MHz_samBA',num2str(ba),...
                '_att0p1f2inc0p',num2str(alpha,'%02d'),'_nc10_80kPa'];
            file_HP = ['PWNE',freq,'MHz_samBA',num2str(ba),...
                '_att0p1f2inc0p',num2str(alpha,'%02d'),'_nc10_400kPa'];


            fnumber = 2;
            file_out = ['RFfn',num2str(fnumber),'_PWNE',freq,...
                'MHz_samBA',num2str(ba),...
                '_att0p1f2inc0p',num2str(alpha,'%02d'),'_nc10_400kPa'];

            disp(file_out)


            tic
            load(fullfile(dataDir,file_LP))
            z = (1:length(rf_prebf))*(1/fs)*c0/2;
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

        end

        clear rf1 rf2

    end
end
%%
dataDir = 'P:\smerino\NonlinAttESt\data\ref_att1p2';

x = linspace(0,0.0384,128); x = x-mean(x);
ba = 10;

for alpha = [40,60]
    for freq = 4:6

        file_LP = ['PWNE',num2str(freq),'MHz_refBA',num2str(ba),...
            '_att0p',num2str(alpha),'f1p2_nc10_80kPa'];
        file_HP = ['PWNE',num2str(freq),'MHz_refBA',num2str(ba),...
            '_att0p',num2str(alpha),'f1p2_nc10_400kPa'];
        fnumber = 2;
        file_out = ['RFfn',num2str(fnumber),'_PWNE',num2str(freq),...
            'MHz_refBA',num2str(ba),...
            '_att0p',num2str(alpha),'f1p2_nc10_400kPa.mat'];
        disp(file_out)


        tic
        load(fullfile(dataDir,file_LP))
        z = (1:length(rf_prebf))*(1/fs)*c0/2;
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

    end

    clear rf1 rf2

end
