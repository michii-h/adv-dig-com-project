clear all
close all
%% Read data from Soundcard
fs = 12000; % sampling rate of base band signal (OFDM)
fs_if = 48000; % Sampling frequency
mTime = 10; % Capturetime in s
r = audiorecorder(fs_if, 16, 1);
recordblocking(r,mTime);
vfCaptureBuffer = getaudiodata(r); % get data
figure(1)
plot(vfCaptureBuffer)
xlabel('\rightarrow k in samples')
ylabel('\rightarrow v(k)')
title('time signal')
%% Evaluate IF Frequency from pilots
% Plot power density spectrum
pwelch(vfCaptureBuffer,[],[],[],fs_if)
% Read from spectrum
f_IF = 1e3*[9.4750 9.4753 9.4750]; % TenTec
f_IF = 1.2001e+04; % DreamTransmitter SET to 9475Hz!!!
%% Digital Down Conversion (DDC)
% plot FFT
V = fft(vfCaptureBuffer);
figure(2)
subplot(4,1,1)
fAxis = (0:length(V)-1)*fs_if/length(V);
plot(fAxis,abs(V))
title('recorded signal (real valued) f_s=48kHz')

% generate analytical signal
V(length(V)/2:end)=0; % set negativ freqnencies to zero
V_analytical = V*2; % same power
subplot(4,1,2)
fAxis = (0:length(V)-1)*fs_if/length(V);
plot(fAxis,abs(V_analytical))
title('analytical signal (complex valued) f_s=48kHz')
v_analytical = ifft(V_analytical);

% shift signal to baseband
%f_IF = mean(f_IF);
f_IF = f_IF(1);
OmegaIF = 2*pi*f_IF/fs_if;
phaser = exp(-1i*OmegaIF*[0:length(v_analytical)-1].');

v_bb = v_analytical.*phaser;

subplot(4,1,3)
fAxis = (0:length(V)-1)*fs_if/length(V);
plot(fAxis,abs(fft(v_bb)))
title('frequency shifted baseband signal (complex valued) f_s=48kHz')

% resample
v_bb = resample(v_bb,1,4);

subplot(4,1,4)
fAxis = (0:length(v_bb)-1)*fs/length(v_bb);
plot(fAxis,abs(fft(v_bb)))
title('resampled baseband signal (complex valued) f_s=12kHz')

%% Detect Robustnes Mode
figure(3)
r_vv = xcorr(v_bb);
vKappa = [0:length(r_vv)-1]-length(r_vv)/2+0.5;

plot(vKappa,abs(r_vv))

xlabel('\rightarrow \kappa')
ylabel('\rightarrow r_{vv}(\kappa)')

RobMetric = [circshift(v_bb,288,1)' ;...
             circshift(v_bb,256,1)' ;...
             circshift(v_bb,176,1)' ;...
             circshift(v_bb,112,1)' ] * v_bb;
         
[~, Mode] = max(abs(RobMetric));


