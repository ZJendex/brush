function [P1, P_plot, f, frange, peakInHz, peakIndex] = plot_powerSpectrum(signal, fs, if_plot)
if nargin < 3 || isempty(if_plot)
    if_plot = 0;
end

% Y = fft(signal);
% L = length(signal);

L = length(signal);
signal = reshape(signal,[1,L]);
w = hamming(L);
Y = fft(signal'.*w);
P2 = abs(Y/L);
P1 = P2(1:floor(L/2)+1);
P1(2:end-1) = 2*P1(2:end-1);

f = fs*(0:floor(L/2))/L;
% figure(17)
frange = (f>0.2) & (f<10000);
P_plot = P1(frange);
% P_plot = P_plot / max(P_plot);


% get the highest peak
% [mV, index] = max(P_plot);

[peaks, locs] = findpeaks(P_plot);
if length(peaks) >= 2
    index = locs(2);
elseif length(peaks) >= 1
    index = locs(1);
else
    index = 0;
end
[vv,index_start_from] = max(frange);
peakInHz = f(index_start_from+index-1);
peakIndex = index_start_from+index-1;

if if_plot
% figure
plot(f(frange),P_plot); hold on
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')
% figure
% plot(f,P1); hold on
% title('Single-Sided Amplitude Spectrum of X(t)')
% xlabel('f (Hz)')
% ylabel('|P1(f)|')
% figure
% plot((1:length(Y_i))/fs, unwrap(angle(Y_i)));
% xlabel('time in second')
% ylabel('phase in rad')
end
end
