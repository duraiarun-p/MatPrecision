clc;
clear;
%%
filename_xlsx = 'Closed_loop.xlsx';
fid_xlsx = fopen(filename_xlsx,'r+');
[a1,b1] = xlsread(filename_xlsx,'Sheet1');
[c1,d1] = xlsread(filename_xlsx,'Sheet2');

Time1 = a1(:,1); Accel11 = a1(:,2)*10; Accel12 = a1(:,3)*10; Accel13 = a1(:,4)*10; Accel14 = a1(:,5)*10;    % Accel coresponding to Mode-1
Time2 = c1(:,1); Accel21 = c1(:,2)*10; Accel22 = c1(:,3)*10; Accel23 = c1(:,4)*10; Accel24 = c1(:,5)*10;    % Accel coresponding to Mode-2

Accel11_part1 = Accel11(1:50000,:);
Accel11_part2 = Accel11(80000:end,:);

%% To plot PSD from time response plot corresponding to Mode-1
Fs = 200; % Sampling frequency 
T = 1/Fs; % Sample time 
L1 = 191522; % Length of signal 
t1 = (0:L1-1)*T; % Time vector
x1 = Accel11;       % Acceleration signals (Mode1)
xdft1 = fft(x1);

figure();
Pxx1 = 1/(L1*Fs)*abs(xdft1(1:length(x1)/2+1)).^2;
freq = 0:Fs/L1:Fs/2;
plot(freq,10*log10(Pxx1));
xlabel('Frequency in Hz'); ylabel('Power/frequency (dB/Hz)');
title('Mode 1')
%% 

%% To plot PSD from time response plot corresponding to Mode-2
Fs = 200; % Sampling frequency 
T = 1/Fs; % Sample time 
L2 = 182834; % Length of signal
t2 = (0:L2-1)*T; % Time vector
x2 = Accel21;       % Acceleration signals (Mode2)
xdft2 = fft(x2);

figure();
Pxx2 = 1/(L2*Fs)*abs(xdft2(1:length(x2)/2+1)).^2;
freq = 0:Fs/L2:Fs/2;
plot(freq,10*log10(Pxx2));
xlabel('Frequency in Hz'); ylabel('Power/frequency (dB/Hz)');
title('Mode 2')
%% 

figure(); plot(Time1,Accel11);
title('Active Vibration Control - Mode-1')
xlabel('Time in Seconds'); ylabel('Right HT tip Aceleration (g)'); grid on
figure(); plot(Time1,Accel13);
title('Active Vibration Control - Mode-1')
xlabel('Time in Seconds'); ylabel('Left HT tip Aceleration (g)'); grid on

figure(); plot(Time2,Accel21);
title('Active Vibration Control - Mode-2')
xlabel('Time in Seconds'); ylabel('Right HT tip Aceleration (g)'); grid on
figure(); plot(Time2,Accel23);
title('Active Vibration Control - Mode-2')
xlabel('Time in Seconds'); ylabel('Left HT tip Aceleration (g)'); grid on
pxx = pcov(Accel11);

% % % % %% Plotting PSD for all the cases (Control On and Control Off)
% % % % Fs = 0.0001;
% % % % % t = 0:1/Fs:1-1/Fs;
% % % % t = 0:(1/Fs):100;
% % % % % x = cos(2*pi*100*t)+randn(size(t));
% % % % % plot(psd(spectrum.periodogram,x,'Fs',Fs,'NFFT',length(x)));
% % % % figure;
% % % % plot(psd(spectrum.periodogram,Accel11,'Fs',Fs,'NFFT',length(Accel11)));
% % % % title('Full Periodogram')
% % % % figure;
% % % % plot(psd(spectrum.welch,Accel11,'Fs',Fs,'NFFT',length(Accel11)));
% % % % title('Fll Welch PSD')
% % % % 
% % % % figure;
% % % % plot(psd(spectrum.periodogram,Accel11_part1,'Fs',Fs,'NFFT',length(Accel11_part1)));
% % % % title('Periodogram with Control ON')
% % % % figure;
% % % % plot(psd(spectrum.welch,Accel11_part1,'Fs',Fs,'NFFT',length(Accel11_part1)));
% % % % title('Welch PSD with Control ON')
% % % % 
% % % % figure;
% % % % plot(psd(spectrum.periodogram,Accel11_part2,'Fs',Fs,'NFFT',length(Accel11_part2)));
% % % % title('Periodogram with Control Off')
% % % % figure;
% % % % plot(psd(spectrum.welch,Accel11_part2,'Fs',Fs,'NFFT',length(Accel11_part2)));
% % % % title('Welch PSD with Control Off')
% % % % %% 

