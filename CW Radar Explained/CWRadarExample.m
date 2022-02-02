% Bryan Tsang

clear;
clc;
close all;

% Read Wav File
info = audioinfo('C:\Users\bryan\OneDrive - The Pennsylvania State University\Desktop\CW Radar Explained\drone.wav');
[rxsignal,fs] = audioread('C:\Users\bryan\OneDrive - The Pennsylvania State University\Desktop\CW Radar Explained\drone.wav');


wv = 0.03;                              % WaveLength (c/freq)
fsamp = 8e3;                            % Desired Sampling Rate
dt = 1/fsamp;                           % Time between samples
Z = 8001;                                  % Starting Index
size = 0.2*fsamp;                         % Specified Integration Time (Total number of Samples = seconds*sampling rate)

rxsignala = resample(rxsignal,fsamp,fs);     % Resample from Audacity's 44.1 KHz to fsamp specified
rxsignalin = rxsignala(:,1)';                % Total Inphase signal
rxsignalq = rxsignala(:,2)';                 % Total Quadrature signal

rxsegin = rxsignalin(Z+1: Z+size+1);        % Inphase signal w/ Specified Intergration Time
rxsegq = rxsignalq(Z+1: Z+size+1);          % Quadrature signal w/ Specified Integration Time
rxseg = (rxsegin + 1i*rxsegq);         % Time Signal = (Inphase + i*Quadrature)


%% Time Signal Plot

N = length(rxseg);
t = 0:dt:(N*dt)-dt;                         % Time axis given # of samples

figure(1);
subplot(2,1,1);
plot(t,rxsegin, 'r')
title('Inphase Time Signal')
xlabel('Time [s]')

subplot(2,1,2);
plot(t,rxsegq, 'b')
title('Quadrature Time Signal')
xlabel('Time [s]')


%% FFT Plot

Y = fft(ifftshift(rxseg))/N;                        % Perform FFT
ff = [0 1:N/2 -N/2:-1] * fsamp/N;                   % Frequency Axis

Y(1:4) = 0;                                         % Filtering out the overpowering low doppler signals
Y(size - 4: size+1) = 0;

figure(2);
plot(ff,abs(Y))
title('10 GHz Power Spectral Density')
xlabel('Freq (Hz)')



%% Spectrogram 
NFFT = 90;                                                      % Spectrogram Specifics
f_vec = [-floor(NFFT/2) : ceil(NFFT/2)-1] * fsamp/NFFT;

figure(3);
% [s,f2,t] = spectrogram(rxseg,NFFT,80,f_vec,fsamp,'yaxis');      % Original
% imagesc(1e3*t,f2*(wv/2),10*log10(abs(s)))                       % Original

[s,f2,t] = spectrogram(rxseg,NFFT,80,'centered');                  % Normalized Freqs
imagesc(t,f2/(pi),10*log10(abs(s)))                                   % Normalized

ax = gca;
ax.YDir = 'normal';      
set(gca,'Visible','off')    %Removes everything except the plot
% z = colorbar;
% title('Velocity vs Time Spectrogram')
% xlabel('Time [ms]')
% % ylabel('Velocity [m/s]')            %Original
% 
% ylabel('Normalized Velocity')                                         %Normalized
% ylabel(z,'Return in dB')
