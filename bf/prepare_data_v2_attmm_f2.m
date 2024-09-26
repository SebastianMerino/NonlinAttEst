clear; close all; clc;
dataDir = 'C:\Users\sebas\Documents\Data\Nonlinearity\attInc\BAinc_ATTinc_4_6_MHz';

fs = 8.555555555555555e+07;
c0 = 1540;

x = linspace(0,0.0384,128); x = x-mean(x);

for freq = 6:6
for attref = 8:2:22
%%
file_LP = ['PWNE',num2str(freq),'MHz_samincBA6inc12_att0p1f2inc0p',num2str(attref,'%02d'),'_nc10_80kPa'];
file_HP = ['PWNE',num2str(freq),'MHz_samincBA6inc12_att0p1f2inc0p',num2str(attref,'%02d'),'_nc10_400kPa'];
% file_LP = 'PWNE4MHz_samBA12_att0p1f2inc0p08_nc10_80kPa';
% file_HP = 'PWNE4MHz_samBA12_att0p1f2inc0p08_nc10_400kPa';
disp(file_LP)
fnumber = 2;

file_out = ['RFfn',num2str(fnumber),'_',file_HP];
%file_out = ['rfLINEAR_',file_HP];

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
% for basam = 6:3:12
% 
%     disp(basam)
% file_LP = ['PWNE_samBA',num2str(basam),'_att0p10f2_nc10_80kPa'];
% file_HP = ['PWNE_samBA',num2str(basam),'_att0p10f2_nc10_400kPa'];
% 
% fnumber = 3;
% 
% file_out = ['RFfn',num2str(fnumber),'_',file_HP];
% %file_out = ['rfLINEAR_',file_HP];
% 
% %keyboard
% load(file_LP)
% z = (1:length(rf_prebf))*(1/fs)*c0/2;
% for rr = 1:size(rf_prebf,3)
%     rf1(:,:,rr) = bf_planewave((rf_prebf(:,:,rr))',fs,fnumber);
% end
% 
% %keyboard
% %rf1 = rf;
% load(file_HP)
% %rf2 = rf;
% for rr = 1:size(rf_prebf,3)
%     rf2(:,:,rr) = bf_planewave((rf_prebf(:,:,rr))',fs,fnumber);
% end
% 
% %keyboard
% save(file_out,'c0','fs','rf1','rf2','x','z');
% 
% end
% 
% clear rf1 rf2
