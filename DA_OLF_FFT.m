clc;
clear;
%%
% dataset = xlsread('60mps_OLF.xlsx');
% dataset = xlsread('Format_All_Original.xlsx','60mps');
% dataset1 = xlsread('Closed_loop.xlsx','Sheet1');
% dataset2 = xlsread('Closed_loop.xlsx','Sheet2');
% 
% Time1 =dataset1(1:end,1);
% Time2 =dataset2(1:end,1);
%%
% Acceleration1 = (dataset1(1:end,2)-mean(dataset1(1:end,2)))/std(dataset1(1:end,2));
% Acceleration2 = (dataset2(1:end,2)-mean(dataset2(1:end,2)))/std(dataset2(1:end,2));
% Fs=roundn(1/(Time1(12,1)-Time1(11,1)),1);
% % Acceleration2 = dataset2(1:end,2);
% % Displacement1 = dataset(2:end,2);
% %%
% wind=50:50:400;
% for i=1:length(wind)
% Acceleration1s=smooth(Acceleration1,wind(i));
% AccelF=fft(Acceleration1s,Fs);
% AccelFa=abs(AccelF);
% figure(1),
% subplot(131),
% plot(AccelFa(1:100)); hold on;
% subplot(132),
% plot(Acceleration1(1:1000)); hold on;
% subplot(133),
% plot(Acceleration1s(1:1000)); hold on;
% pause(0.5);
% end
%%
filename_xlsx = 'Closed_loop.xlsx';
fid_xlsx = fopen(filename_xlsx,'r+');
[a1,b1] = xlsread(filename_xlsx,'Sheet1');
[c1,d1] = xlsread(filename_xlsx,'Sheet2');

Time1 = a1(:,1); Accel11 = a1(:,2)*10; Accel12 = a1(:,3)*10; Accel13 = a1(:,4)*10; Accel14 = a1(:,5)*10;    % Accel coresponding to Mode-1
Time2 = c1(:,1); Accel21 = c1(:,2)*10; Accel22 = c1(:,3)*10; Accel23 = c1(:,4)*10; Accel24 = c1(:,5)*10;    % Accel coresponding to Mode-2
%%
Fs1=roundn(1/(Time1(12,1)-Time1(11,1)),1);
Fs2=roundn(1/(Time2(12,1)-Time2(11,1)),1);
%%
smoothwindowN=300;
lastfreqofoneside=200;
[A11,f1]=onesideFFT(Accel11,Fs1,smoothwindowN,lastfreqofoneside);
[A12,~]=onesideFFT(Accel12,Fs1,smoothwindowN,lastfreqofoneside);
[A13,~]=onesideFFT(Accel13,Fs1,smoothwindowN,lastfreqofoneside);
[A14,~]=onesideFFT(Accel14,Fs1,smoothwindowN,lastfreqofoneside);

[A21,f2]=onesideFFT(Accel21,Fs2,smoothwindowN,lastfreqofoneside);
[A22,~]=onesideFFT(Accel22,Fs2,smoothwindowN,lastfreqofoneside);
[A23,~]=onesideFFT(Accel23,Fs2,smoothwindowN,lastfreqofoneside);
[A24,~]=onesideFFT(Accel24,Fs2,smoothwindowN,lastfreqofoneside);
%% signal plot
figure(1),
subplot(2,2,1),
plot(f1,A11);
subplot(2,2,2),
plot(f1,A12);
subplot(2,2,3),
plot(f1,A13);
subplot(2,2,4),
plot(f1,A14);

figure(2),
subplot(2,2,1),
plot(f2,A21);
subplot(2,2,2),
plot(f2,A22);
subplot(2,2,3),
plot(f2,A23);
subplot(2,2,4),
plot(f2,A24);
%% Function one-sided FFT

function [P1,f]=onesideFFT(X,Fs,smoothwindowN,lastfreqofoneside)
X=(X-mean(X))/std(X); % Zero mean - Zero variance (Normalisation) - Zero baseline
X=smooth(X,smoothwindowN); % Smoothening data with smooth function that is a low-pass filter, choose window size based on sampling frequency
L=length(X);
Y = fft(X); % FFT at full signal length
P2 = abs(Y/L); % magnitude 
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1); % one-sided power spectrum
f = Fs*(0:(L/2))/L; % frequency bin generated upto Fs/2
[val,idx]=min(abs(f-lastfreqofoneside)); % User defined frequency for one-sided spectrum
f=f(1:idx);
P1=P1(1:idx);
end