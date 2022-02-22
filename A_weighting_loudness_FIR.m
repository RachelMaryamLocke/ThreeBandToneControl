%FIR filter design frequency - response method for loudness control
%(inverse A - weighted curve)

close all
clear all

Fs = 44100; %sampling frequency 
N = 2048;%DFT length
f = (0:N/2-1)'/(N/2) * Fs/2; %Creates a frequency vector

%%
%A - weighting curve equation

fSq = f .* f;

H = ((fSq + 20.6 .^ 2) .* (fSq + 12194 .^ 2) .* sqrt((fSq + 107.7 .^ 2).*(fSq + 737.9 .^ 2))) ./ ((12194 .^ 2) .* f .^ 4);
H(1) = 0; 

%tests A-weighting curve
figure (1);
subplot(211)
semilogx(f, (20*log10(H)),'linewidth', 2);
xlabel('Frequency [Hz]');
ylabel('Amplitude [dB]');
title('Inverse A-Weighting Curve');
ylim([0 52])
xlim([0 30000])

% %%
% % Symmetry
% 
% H = [H; flipud(H(2:end))]; %alternative (to below) for H using coloums instead of rows 
% % (uses one line instead of two and avoids concatination errors)
%  
% % R(1:1024) = InvA;
% % R(22051:44100 - 1) = fliplr(InvA(2:end));
% 
% %tests symmetry
% subplot(212)
% plot(20*log10(H),'linewidth',2) %xaxis frequency Hz (logarithimic), yaxis dB (20*log10)
% grid on
% xlabel('samples[N]')
% ylabel('Gain[dB]')
% title('A-weighting Curve Symmetry')
% xlim([0 2100])
% %%
% 
% %take the inverse fast fourier transform to obtain h(n) 
%  
% h = ifft(H); %inverse FFT of symmetrical impulse response
% 
% 
% figure(2)
% subplot(211)
% stem(abs(real(h)));
% title('Impulse Response Before fftshift')
% ylabel('Amplitude')
% xlabel('Time')
% xlim([0 2000])
% 
% h = fftshift(h); %make it causal
% 
% %plots causal impulse response
% subplot(212)
% stem(abs(h),'g-o')
% title('Impulse Response After fftshift')
% ylabel('Amplitude')
% xlabel('Time')
% xlim([0 2000])
% 
% %%
% %time limit the impulse responce
% 
% M = 100;
% trunc = round(length(h)/2) + (-M:M);
% ht = h(trunc);
% 
% figure(3)
% subplot(211)
% stem(abs(h),'g-o')
% grid on
% title('Impulse Response After fftshift')
% ylabel('h[n]')
% xlabel('time[s]')
% xlim([0 2000])
% 
% subplot(212)
% plot(1:length(h), h, trunc, ht, 'r-o')
% grid on
% title('Impulse Response After Truncation Before Windowing')
% xlabel('time[s]')
% ylabel('h[n]')
% ylim([-0.5 3])
% xlim([900 1150])
% 
% %%
% % window the impulse response
% 
% h = ht.*hann(length(ht)); %window the signal to obtain the windowed impulse response 
% 
% figure(4)
% subplot(211)
% plot(ht,'o-r')
% grid on
% title('Truncated Impulse Response Before Windowing')
% xlabel('time[s]')
% ylabel('h[n]')
% xlim([0 210])
% 
% subplot(212)
% plot(h,'o-b')
% grid on
% title('Truncated Impulse Response After Windowing')
% xlabel('time[s]')
% ylabel('h[n]')
% xlim([0 210]) 
% %%
% 
% %FIR filter coefficients 
% 
% %The values of the impulse response of a non-recursive
% %filter are the filter coefficients themselves
% 
% bn = h;
% 
% %%
% 
% %test the filter
% 
% % [X, Fs] = audioread('longer-sample');
% 
% x = randn(10000, 1) * 2 - 1; %create a random signal to test the filter
% xf = filter(bn, 1, x); %filter the signal 
% 
% % x = zeros(10000, 1);%check it's working by feeding an impulse to the filter
% % x(1) = 1;
% 
% F = (0:length(x)-1)/length(x)*Fs;%create a freqency axis
% figure(5)
% semilogx(F, 20*log10(abs(fft([x xf])))) %plots frequency axis F against, the FFT of the original signal and the filtered signal
% grid on
% title('FFT of Noise Signal and Filtered Noise')
% legend('Unfiltered Noise', 'Filtered Noise')
% ylabel('Amplitude[dB]')
% xlabel('Frequency[Hz]')
% xlim([20 20000])
% 
% 
% %%