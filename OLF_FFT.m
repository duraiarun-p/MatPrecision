clc
clear

% dataset = xlsread('60mps_OLF.xlsx');
% dataset = xlsread('Format_All_Original.xlsx','60mps');
dataset1 = xlsread('Closed_loop.xlsx','Sheet1');
dataset2 = xlsread('Closed_loop.xlsx','Sheet2');

Time1 =dataset1(1:end,1);
Time2 =dataset2(1:end,1);
Acceleration1 = dataset1(1:end,2);
Acceleration2 = dataset2(1:end,2);
% Displacement1 = dataset(2:end,2);
% Displacement2 = dataset(2:end,4);

%% FFT from acceleration

fs = 1/0.001; % sampling frequency
LL1 = length(Time1);
LL2 = length(Time2);
N=4096;  % Block size
f = fs/2*linspace(0, 1, N/2+1);  % Frequency values
%%
% subplot(211)
figure()
X = fft(Acceleration1, N)/LL1;          % Accel1 on Right HT
plot(f, 2*abs(X(1:N/2+1)));  grid on; 
xlabel('Frequency in Hz');
ylabel('Amplitude');
title('One-sided Amplitude Spectrum') , shg
%%
% subplot(212)
figure()
Y = fft(Acceleration2, N)/LL2;          % Accel2 on Left HT   
plot(f, 2*abs(Y(1:N/2+1)));  grid on; 
xlabel('Frequency in Hz');
ylabel('Amplitude');
%%
figure()
plot(f, 2*abs(X(1:N/2+1)), 'r', f, 2*abs(Y(1:N/2+1)), 'g');  grid on; 
xlabel('f, frequency, [Hz]')
ylabel('Amplitude')
legend('Acceleration Response - Mode 1 ','Acceleration Response - Mode 2')
title('One-sided Amplitude Spectrum') , shg
%%
% figure(3)
% plot(Time1,Displacement1, 'r', Time1, Displacement2, 'g');

figure()
plot(Time1,Acceleration1);
xlabel('Time in seconds')
ylabel('Acceleration in g')
legend('Time Response-Mode1')
%%
figure()
plot(Time2,Acceleration2);
xlabel('Time in seconds')
ylabel('Acceleration in g')
legend('Time Response-Mode2')
%%
% % Fs = 0.0001;
% % % t = 0:1/Fs:1-1/Fs;
% % t = 0:(1/Fs):100;
% x = cos(2*pi*100*t)+randn(size(t));
% % plot(psd(spectrum.periodogram,x,'Fs',Fs,'NFFT',length(x)));
% figure;
% plot(psd(spectrum.periodogram,Acceleration1,'fs',fs,'NFFT',length(Acceleration1)));
% title('Full Periodogram - Mode1')
% figure;
% plot(psd(spectrum.welch,Acceleration1,'fs',fs,'NFFT',length(Acceleration1)));
% title('Full Welch PSD - Mode1')
% 
% figure;
% plot(psd(spectrum.periodogram,Acceleration2,'fs',fs,'NFFT',length(Acceleration2)));
% title('Full Periodogram - Mode2')
% figure;
% plot(psd(spectrum.welch,Acceleration2,'fs',fs,'NFFT',length(Acceleration2)));
% title('Full Welch PSD - Mode2')

