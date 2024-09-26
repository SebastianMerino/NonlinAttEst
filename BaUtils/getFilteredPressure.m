function  P = getFilteredPressure(rf,fs,freqC,freqTol,order)
%   Filters the rf data around a center frequency and returns the envelope
%   All frequencies should be in MHz

t = (0:size(rf,1)-1)'./fs;
rfMod = rf.*exp(1j*2*pi*freqC*t);

% if plotS
% Srf = db(mean(abs(fft(rf)),2));
% Srf = Srf - max(Srf);
% f = (0:size(rf,1)-1)'./size(rf,1)*fs;
% figure, plot(f,Srf)
% xline(freqC-freqTol)
% xline(freqC+freqTol)
% xlim([0 2*freqC])
% end

freqNyq = fs/2;
d = floor(order/2);
b = fir1(order,freqTol/freqNyq);
[~,n,p] = size(rf);
rfFilt = filter(b,1,[rfMod;zeros(d,n,p)]);
rfFilt = rfFilt(d+1:end,:,:);
P = abs(rfFilt);

% Srf = fftshift(db(mean(abs(fft(rfFilt)),2)));
% Srf = Srf - max(Srf);
% f = (0:size(rf,1)-1)'./size(rf,1)*fs - fs/2;
% figure, plot(f,Srf)
% xline(-freqTol)
% xline(freqTol)
% xlim([-freqC freqC])

end