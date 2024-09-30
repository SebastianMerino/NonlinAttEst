clear; close all; clc;
dataDir = 'C:\Users\sebas\Documents\Data\Nonlinearity\homoAttMismatch';

fs = 8.555555555555555e+07;
c0 = 1540;

x = linspace(0,0.0384,128); x = x-mean(x);

for freq = 5:5
for attref = 8:2:22
%%
% file_LP = ['PWNE',num2str(freq),'MHz_refBA6_att0p',num2str(attref,'%02d'),'f2p0_80kPa'];
% file_HP = ['PWNE',num2str(freq),'MHz_refBA6_att0p',num2str(attref,'%02d'),'f2p0_400kPa'];
% file_out = ['RFfn',num2str(fnumber),'_',file_HP];

file_LP = ['PWNE_refBA6_att',num2str(attref),'f2_nc10_80kPa'];
file_HP = ['PWNE_refBA6_att',num2str(attref),'f2_nc10_400kPa'];
fnumber = 2;
file_out = ['RFfn',num2str(fnumber),'_PWNE',num2str(freq),'MHz_refBA6_att0p',num2str(attref,'%02d'),'f2p0_400kPa'];
disp(file_out)


%keyboard
load(fullfile(dataDir,file_LP))
z = (1:length(rf_prebf))*(1/fs)*c0/2;
for rr = 1:size(rf_prebf,3)
    rf1(:,:,rr) = bf_planewave((rf_prebf(:,:,rr))',fs,fnumber);
end

%keyboard
%rf1 = rf;
load(fullfile(dataDir,file_HP))
%rf2 = rf;
for rr = 1:size(rf_prebf,3)
    rf2(:,:,rr) = bf_planewave((rf_prebf(:,:,rr))',fs,fnumber);
end

%keyboard
save(fullfile(dataDir,file_out),'c0','fs','rf1','rf2','x','z');

end

clear rf1 rf2

end

%%
% for basam = 6:3:12
% 
% %    disp(basam)
% file_LP = ['PWNE4MHz_samBA',num2str(basam),'_att0p1f2_nc10_80kPa'];
% file_HP = ['PWNE4MHz_samBA',num2str(basam),'_att0p1f2_nc10_400kPa'];
% fnumber = 2;
% file_out = ['RFfn',num2str(fnumber),'_PWNE4MHz_samBA',num2str(basam),'_att0p10f2_nc10_400kPa'];
% disp(file_out)
% 
% %keyboard
% load(fullfile(dataDir,file_LP))
% z = (1:length(rf_prebf))*(1/fs)*c0/2;
% for rr = 1:size(rf_prebf,3)
%     rf1(:,:,rr) = bf_planewave((rf_prebf(:,:,rr))',fs,fnumber);
% end
% 
% %keyboard
% %rf1 = rf;
% load(fullfile(dataDir,file_HP))
% %rf2 = rf;
% for rr = 1:size(rf_prebf,3)
%     rf2(:,:,rr) = bf_planewave((rf_prebf(:,:,rr))',fs,fnumber);
% end
% 
% %keyboard
% save(fullfile(dataDir,file_out),'c0','fs','rf1','rf2','x','z');
% 
% end
% 
% clear rf1 rf2
