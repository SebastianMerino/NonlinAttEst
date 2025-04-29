%% Data Block test
% How much does data varies within a data block?
clear, clc
gamma = 1.4;
f = (3:0.1:7)'; % 3 to 7 MHz

%% Variation per depth
z = 0.5:0.5:5;

bona = 9;
b = 1 + bona/2;
a = 0.15/db(exp(1)); % Attenuation in Np/cm/MHz^gamma

data = b.*(1-exp(-2*a.*f.^gamma.*z))./(a.*f.^gamma.*z);
figure,
plot(f,data)

%% Variation per nonlinearity
bona = 6:12;
b = 1 + bona/2;
a = 0.15/db(exp(1)); % Attenuation in Np/cm/MHz^gamma
z = 3; % depth in cm

data = b.*(1-exp(-2*a.*f.^gamma.*z))./(a.*f.^gamma.*z);
figure,
plot(f,data)

%% Variation per attenuation
bona = 9;
b = 1 + bona/2;
a = (0.08:0.02:0.16)/db(exp(1)); % Attenuation in Np/cm/MHz^gamma
z = 3; % depth in cm

data = b.*(1-exp(-2*a.*f.^gamma.*z))./(a.*f.^gamma.*z);
figure,
plot(f,data)

%% Variating both: same ratio
goldratio = (1+9/2)/0.05;
a = (0.08:0.02:0.16)/db(exp(1)); % Attenuation in Np/cm/MHz^gamma
b = goldratio*a;
z = 3;

data = b.*(1-exp(-2*a.*f.^gamma.*z))./(a.*f.^gamma.*z);
figure,
plot(f,data)

%% ==================================================================== %%
%% Stability of system
clear, clc
gamma = 2;
f(1,1,:) = (3:0.5:7)';

bona0 = 9;
b0 = 1 + bona0/2;
a0 = 0.14/db(exp(1)); % Attenuation in Np/cm/MHz^gamma
z = 3;

for i=1:5
data = b0.*(1-exp(-2*a0.*f.^gamma.*z))./(a0.*f.^gamma.*z);
data = data + 0.1*randn(size(data));

a = (0.08:0.001:0.16)'/db(exp(1)); % Attenuation in Np/cm/MHz^gamma
bona = 6:0.05:11;
b = 1 + bona/2;

model = b.*(1-exp(-2*a.*f.^gamma.*z))./(a.*f.^gamma.*z);
l2norm = sum((model-data).^2,3)/2;
figure("Position",[100 100 300 250]),
imagesc(bona,a*db(exp(1)),db(l2norm))
ax = gca;
ax.YDir = "normal";
colorbar
xlabel('B/A')
ylabel('AC [dB/cm/MHz^\gamma]')
title("L2 norm")
hold on
goldberg = b0/a0;
aRatio = b/goldberg;
% plot(bona,aRatio*db(exp(1)), "w--", 'LineWidth',2)
plot(bona0,a0*db(exp(1)), "wx", 'LineWidth',2)
hold off

end
%% Stability of previous system
clear, clc
gamma = 2;
f = 5;
z(1,1,:) = 3:0.05:3.5;

bona0 = 9;
b0 = 1 + bona0/2;
a0 = 0.14/db(exp(1)); % Attenuation in Np/cm/MHz^gamma

for ii=1:5
data = b0.*(1-exp(-2*a0.*f.^gamma.*z))./(a0.*f.^gamma.*z);
data = data + 0.1*randn(size(data));

a = (0.08:0.001:0.16)'/db(exp(1)); % Attenuation in Np/cm/MHz^gamma
bona = 6:0.05:11;
b = 1 + bona/2;

 tg-ñl = b.*(1-exp(-2*a.*f.^gamma.*z))./(a.*f.^gamma.*z);
l2norm = sum((model-data).^2,3)/2;

figure("Position",[100 100 300 250]),
imagesc(bona,a*db(exp(1)),db(l2norm))
ax = gca;
ax.YDir = "normal";
colorbar
xlabel('B/A')
ylabel('AC [dB/cm/MHz^\gamma]')
title("L2 norm")
hold on
goldberg = b0/a0;
aRatio = b/goldberg;
plot(bona,aRatio*db(exp(1)), "w--", 'LineWidth',2)
plot(bona0,a0*db(exp(1)), "wx", 'LineWidth',2)
hold off
end
