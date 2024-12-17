%% Capture Audio Signal
% Generate Object
fs = 48000; % Sampling frequency of if signal
recObj = audiorecorder(fs,16,1);
recordblocking(recObj, 5);
vfCaptureBuffer = getaudiodata(recObj);
figure(1)
subplot(2,1,1)
plot(vfCaptureBuffer);
subplot(2,1,2)
pwelch(vfCaptureBuffer,[],[],[],fs);

f_if = 9.4783e+03+0;

%% Digital Down Conversion (ddc)
figure(2)
% analytical signal
S = fft(vfCaptureBuffer);
subplot(4,1,1)
plot(abs(S))

S(length(S)/2:end) = 0;
s_if = ifft(S);
subplot(4,1,2)
plot(abs(S))

% shift from if to baseband
Omega_IF = 2*pi*f_if/fs;
k = [0:length(s_if)-1].';
x = s_if.*exp(-j*Omega_IF*k);
subplot(4,1,3)
plot(abs(fft(x)))

% resample signal from 48kHz to 12kHz
x = resample(x,1,4);
fs = 12000;
subplot(4,1,4)
plot(abs(fft(x)))

figure(3)
pwelch(x,[],[],[],fs)

%% Estimate Robustness Mode
mcX = [circshift(x,[288 0])';...
       circshift(x,[256 0])';...
       circshift(x,[176 0])';...
       circshift(x,[112 0])'];
   
rbMetric = mcX *x;
[value iMode] = max(abs(rbMetric));
