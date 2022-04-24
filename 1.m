%% ELEG 3123 Project - ECG Processing
%  Group:
%  Members:
clear; close all; clc;
% Some parameters for this application
Mfilter_Length = 601;               % Median filter length to estimate baseline
nstart = 101000;                    % Index for beginning of 2 cycle epoch
nstop = 102749;                     % End of 2 cycle epoch
% Select one of nine leads to analysis [1 - 9] = {V1,V2,V3,V4,V5,V6,I,II,III}
Lead=7;
% Data stored in a structure
% Change the address to where you stored the file on your computer
x=load('ECG_Recording1.mat');
% Display the upper level fieldnames
xNames = fieldnames(x);
disp('Structure Names')
disp(xNames)
% Display the fieldnames of Data
xDataNames = fieldnames(x.Data);
disp('Names of Data Structure')
disp(xDataNames)
Fs = x.Data.Sampling_frequency; Fsd2 = Fs/2; % Fs is original sampling freq
Lead_Label = x.Data.Labels{Lead};
%% Extract ECG signal (x.Data.ECG is N=6*60*1000 rows 9 column matrix of ECG Data)
y = x.Data.ECG(:,Lead);
tm=0:1/Fs:(length(y)-1)/Fs;
%plot(tm,y);grid;
%xlabel('time (sec)'); ylabel(['Magnitude Lead ',char(Lead_Label),' (V)']);
%% Remove baseline drift
figure()
% Scale data to mV
y = 1e3*y;
%plot(y);grid;hold on
% Median Filter to remove baseline drift
ymed=medfilt1(y,Mfilter_Length);
plot(tm,y,tm,ymed);grid;
xlabel('time (sec)'); ylabel(['Magnitude Lead ',char(Lead_Label),' (mV)']);
title('ECG Signal with Estimated Baseline')
legend('EKG','Baseline')
set(gca,'Fontsize',14)
set(gca,'Fontsize',14)
y = y-ymed;
figure()
plot(tm,y);grid
xlabel('time (sec)')
ylabel(['Magnitude Lead ',char(Lead_Label),' (mV)']);
title('Baseline Drift Removed')
set(gca,'Fontsize',14)
figure()
plot(tm(nstart:nstop),y(nstart:nstop));grid
xlabel('time (sec)')
ylabel(['Magnitude Lead ',char(Lead_Label),' (mV)']);
title('Baseline Drift Removed Contains 50Hz Noise')
xlim([101 103])
set(gca,'Fontsize',14)

%% Remove Power Line Interference (50Hz)
% Enter the required values for W0 and Bw
W0 = 50/(1000/2); % 50 is the notch filter frequency. Fs is the original sample
% rate, and W0 is measured in Pi radians per second. 
Q = 100000;
Bw = W0/Q; % denominator sets the quality factor. A lower
% Q is snappier, while a higher Q is bouncier.

% MATLAB Notch Filter - Use the following if you are using MATLAB
[b,a]=iirnotch(W0,Bw); 
% Filter data to remove 50Hz line noise
yf = filtfilt(b,a,y);
figure()
plot(tm,yf);grid
xlabel('time (sec)')
ylabel(['Magnitude Lead ',char(Lead_Label),' (mV)']);
title('Power Line (50 Hz) Removed')
xlim([101 103])
set(gca,'Fontsize',14)

%% Low Pass Filter with Fstop = 125 Hz
% Design a FIR low pass filter
Fstop = [125 500];
Fpass = [0 105];
Astop = 60;
Apass = 0.2;

N = 126;         % filter order
F = [Fpass Fstop]/(Fs/2);         % filter frequency edges
A = [1 1 .1 .1];         % filter magnitude at the edges
W = [Apass Astop];         % filter weights one per band
% Matlab - Coefficients for a Lowpass FIR Filter
b = firpm(N,F,A,W);
% Both Matlab and Octave - Coefficients for a Lowpass FIR Filter 
% b = remez(N,F,A,W);   % FIR impulse response (numerator coefficients)
a = 1;                % Denominator multiplier
% Display the Magnitude Frequency Response of the Lowpass FIR Filter
[h,w] = freqz(b,a);   % Filter frequency response
figure()
plot(Fsd2*w/pi,mag2db(abs(h)));grid
xlabel('Frequency (Hz)');ylabel('Magnitude (dB)')
title('Low Pass Equiripple Filter')
set(gca,'Fontsize',14)
% Filter the signal
yf = filtfilt(b,a,yf);
figure()
plot(tm(nstart:nstop),yf(nstart:nstop));grid
xlabel('time (sec)')
ylabel(['Magnitude Lead ',char(Lead_Label),' (mV)']);
title('Low Pass Filtered Signal')
xlim([101 103])
set(gca,'Fontsize',14)

%% Power Spectral Density
L=Fs*10;                % Define segment length
Nfft = 2048;
window = blackmanharris(L,'symmetric');
overlap = L/2;
% Enter the parameters for pwelch to find the PSD of y and yf
[Pyy,f] = pwelch(y,window,overlap, Nfft);
[Pyfyf,f] = pwelch(yf,window,overlap, Nfft);
pf=plot(f,mag2db(Pyy),f,mag2db(Pyfyf));grid
xlabel('Frequency (Hz)');
ylabel('S_{yy}(f) (dB)')
title('Power Spectral Density for Original and Filtered Signals')
legend('Original','Filtered');
set(pf,'LineWidth',2)
set(gca,'Fontsize',14)

%% Reduce the sampling rate of the signal from Fs=1000Hz to FS=250Hz
% enter the code to perform sample rate reduction and display the signal
disp("Reducing");
z1 = decimate(y, 4);
zf1 = decimate(yf,4);
disp("pause")
pause
[Pzz,f] = pwelch(z1, window, overlap, Nfft);
[Pzfzf,f] = pwelch(zf1, window, overlap, Nfft);
figure();
disp("pause")
pause
pf=plot(f,mag2db(Pzz),f,mag2db(Pzfzf));grid
xlabel('Frequency (Hz)');
ylabel('S_{zz}(f) (dB)')
title('Power Spectral Density for Origina and Filtered Signals')
legend('Original','Filtered');
set(pf,'LineWidth',2)
set(gca,'Fontsize',14)
