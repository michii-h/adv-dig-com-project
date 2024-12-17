function [x_baseband] = captureDRM(SwitchCase)

switch SwitchCase
    case 1
        % %% Capture DRM Signal
        fs = 48000;
        recObj = audiorecorder(fs, 16, 2);
        disp('Start speaking.')
        recordblocking(recObj, 5);
        disp('End of Recording.');
        myRecording = getaudiodata(recObj);
    case 2
        % save 'drm.mat' myRecording fs
        load 'drm.mat'
end
% Plot the waveform.
figure(1)
pwelch(myRecording(:,1),[],[],[],48000);

vfCapturebuffer = myRecording(:,1);
figure(2)
subplot(4,1,1)
plot(abs(vfCapturebuffer).^2)


subplot(4,1,2)
L = length(vfCapturebuffer);
f = fs*([0:L-1]-L/2)/L;
plot(f,abs(fftshift(fft(vfCapturebuffer))))

f_IF = 9.4768e+03;

%% Digital Down Conversion (DDC)
% Analytical signal
V = fft(vfCapturebuffer); % Calculate spectrum
V(length(V)/2:end)=0; % set negative frequencies to zero
v_analytical = ifft(2*V);

% Down conversion to baseband by multiplying with phaser
Omega_if = 2*pi*f_IF/fs; % normalized frequency
k = [0:length(v_analytical)-1].'; % sample index
x_baseband = v_analytical .* exp(-j*Omega_if*k); % down conversion

figure(2)
subplot(4,1,3)
L = length(vfCapturebuffer);
f = fs*([0:L-1]-L/2)/L;
plot(f,abs(fftshift(fft(x_baseband))))

% Resampling from 48kHz to 12kHz (required by DRM..mandatory)
x_baseband = resample(x_baseband,1,4); % Downsample by factor 4

figure(2)
subplot(4,1,4)
L = length(x_baseband);
f = 12000*([0:L-1]-L/2)/L;
plot(f,abs(fftshift(fft(x_baseband))))



