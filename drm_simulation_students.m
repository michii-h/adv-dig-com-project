%% OFDM Development Environment
% This simulation implements a complete DRM (Digital Radio Mondiale) system
% with OFDM modulation for image transmission over various channels
start_up;

SwitchDemoSync = false; % Flag to enable additional visualization plots

%% Transmission Parameters Initialization
% The satellite link model uses the following parameters:
% - SNR calculation: $SNR_{dB} = P_{tx} - L_{path} - N_{0} + G_{rx}$
% - Path loss: $L_{path} = 20\log_{10}(4\pi d f/c)$ dB
% - Noise floor: $N_{0} = -174 + 10\log_{10}(BW) + NF$ dBm
tx_power = 7 + 30; % 7 dBm output power of Adalm Pluto + 30 dBm amplifier
link = SatelliteLink('tx_center_frequency', 2400.113 * 1e6, 'tx_gain', 0, 'rx_gain', 30, 'tx_power', tx_power);
link.useSameFrequency = true;
link.displayLinkBudget();

%% Bit Error Rate Calculation from expected SNR
% Theoretical BER for M-QAM in AWGN:
% $P_b \approx \frac{4}{\log_2 M}(1-\frac{1}{\sqrt{M}})\cdot Q(\sqrt{\frac{3\log_2 M}{M-1}\cdot\frac{E_b}{N_0}})$
% where $Q(x) = \frac{1}{\sqrt{2\pi}}\int_{x}^{\infty}e^{-t^2/2}dt$
ModOrder = 16;
BER_sim = berawgn(link.SNR_linear, 'qam', ModOrder);

fprintf('  Simulated BER: %.2e\n', BER_sim);

%% DRM Initialization
% DRM standard defines robustness modes (A-D) and spectrum occupancy (0-5)
% - Mode D: Designed for NVIS (Near Vertical Incidence Skywave) with high robustness
% - Occupancy 3: Bandwidth of ~10 kHz
stDRM.mode = 4;      % Corresponds to Mode D
stDRM.occupancy = 3;

check_drm_bandwidth(link.baseStation.baseband_sample_rate, stDRM.mode, stDRM.occupancy);

%% OFDM Initialization
% OFDM parameters based on DRM standard:
% - $N_{FFT}$: FFT size for the mode/occupancy
% - $N_g$: Guard interval length (cyclic prefix)
% - $N_s = N_{FFT} + N_g$: Total symbol length
stOFDM.iNfft = get_drm_n_useful(stDRM.mode,stDRM.occupancy);
stOFDM.iNg = get_drm_n_guard(stDRM.mode,stDRM.occupancy);
stOFDM.iNs = stOFDM.iNfft + stOFDM.iNg;

%% Generate DRM Frame
% Encodes image data into DRM frames with call sign
% Each frame contains multiple OFDM symbols with data and pilots
image_path = 'th-rosenheim-logo-colored-square.png';
call_sign = 'DL0FHR';
[Slk, M, image_size, iNofFramesNeeded, iNOfFrames] = generate_drm_frames(stDRM, stOFDM, image_path, call_sign);

%% OFDM Modulator
% OFDM modulation process:
% 1. IFFT: $s(n) = \frac{1}{N}\sum_{k=0}^{N-1} S(k) e^{j2\pi kn/N}$
% 2. Add cyclic prefix: copy last $N_g$ samples to beginning
% 3. Parallel to serial conversion
SlkTemp = ifft(fftshift(Slk,2),stOFDM.iNfft,2);

% Add Cyclic Prefix: copy last Ng samples to beginning
% CP prevents ICI and helps with synchronization: $s_{CP}(n) = s(n-N_g)$ for $n=0,1,...,N_g-1$
SlkTemp = [SlkTemp(:,end-stOFDM.iNg+1:end) SlkTemp];

% Parallel to Serial Conversion: reshape to 1-D signal
SlkTemp = SlkTemp.';
vfcTransmitSignal = SlkTemp(:);

%% Channel
% Different channel options:
% 0: Ideal (no distortion)
% 1: Simulated multipath: $y(t) = \sum_{i=0}^{L-1} h_i(t) x(t-\tau_i) + n(t)$
% 2: Simulated QO-100 satellite channel
% 3: Adalm Pluto SDR with satellite correction
% 4: Adalm Pluto SDR loopback
iSwitchChannel = 4;

switch iSwitchChannel

    case 0 % ideal channel
        % Repeat iG Frames
        iG = 2;
        vfcTransmitSignal = repmat(vfcTransmitSignal,iG,1);

        fprintf('Using ideal channel...\n');
        vfcReceiveSignal = vfcTransmitSignal;

    case 1 % simulated channel
        fprintf('Using simulated channel...\n');

        % Repeat iG Frames
        iG = 2;
        vfcTransmitSignal = repmat(vfcTransmitSignal,iG,1);

        stChannel = initChannel();
        vfcReceiveSignal = channel_sim(vfcTransmitSignal,stChannel);

    case 2 % Simulate Satellite Communication for Q0-100
        % Repeat iG Frames
        iG = 2;
        vfcTransmitSignal = repmat(vfcTransmitSignal,iG,1);

        % Update to use the link object instead of stSat
        vfcReceiveSignal = simulate_qo100_channel(vfcTransmitSignal, link);

    case 3 % Use Adalm Pluto
        % SDR come from the link object
        vfcReceiveSignal = LoopbackAdalmPlutoSat(vfcTransmitSignal, link);

        % Frequency offset compensation
        % Manual frequency offset correction: $y_{corr}(n) = y(n) \cdot e^{-j2\pi \Delta f \cdot n/f_s}$
        freqShift = 282; % Measured frequency offset in Hz

        % Convert to angular frequency: $\Delta\Omega = 2\pi\Delta f/f_s$
        delta_Omega = 2*pi*(freqShift/link.baseStation.baseband_sample_rate);
        vPhaser = exp(-1i * delta_Omega * [0:length(vfcReceiveSignal)-1])';

        vfcReceiveSignal = vfcReceiveSignal .* vPhaser;

        % Display spectrum with pilot markers
        figure(375);
        pwelch(vfcReceiveSignal, [],[],[], link.baseStation.baseband_sample_rate);
        hold on;
        xline(link.baseStation.baseband_sample_rate*0.750/12000);
        xline(link.baseStation.baseband_sample_rate*2.250/12000);
        xline(link.baseStation.baseband_sample_rate*3.000/12000);
        xticks([ ...
            link.baseStation.baseband_sample_rate*0.750/12000,...
            link.baseStation.baseband_sample_rate*2.250/12000,...
            link.baseStation.baseband_sample_rate*3.000/12000
            ])

        case 4 % Use Adalm Pluto
        % SDR come from the link object
        vfcReceiveSignal = LoopbackAdalmPluto(vfcTransmitSignal, link);
end

%% Detect Robustness Mode
% Detection algorithm uses autocorrelation between received signal and each possible mode
% $R_{xx}(\Delta k) = \sum_{n=0}^{N-1} x(n) \cdot x^*(n-\Delta k)$
% The mode with highest correlation peak is selected
[iModeEst, iNfft, iNg, iNs, iNOfSymbolsPerFrame] = detect_robustness_mode(vfcReceiveSignal, stDRM.occupancy, SwitchDemoSync);

%% Synchronization
% OFDM synchronization uses cyclic prefix correlation:
% $R(k) = \sum_{n=0}^{N_g-1} x^*(n+k) \cdot x(n+k+N_{FFT})$
% - Time sync: $\hat{k} = \arg\max_k |R(k)|$
% - Frequency offset: $\Delta\hat{f} = \frac{\angle R(\hat{k})}{2\pi T_{useful}}$
vfcReceiveSignal = sync(vfcReceiveSignal, iNs, iNg, iNfft, SwitchDemoSync);

%% OFDM Demodulator
% Reverses the modulation process:
% 1. Serial to parallel conversion
% 2. Remove cyclic prefix
% 3. FFT: $R(k) = \sum_{n=0}^{N-1} r(n) e^{-j2\pi kn/N}$
RlkTemp = reshape(vfcReceiveSignal,stOFDM.iNs,[]).';

% Remove Cyclic Prefix
RlkTemp = RlkTemp(:,stOFDM.iNg+1:end);

% FFT
Rlk = fftshift(fft(RlkTemp,stOFDM.iNfft,2),2);

%% Frame Detection
% Uses cross-correlation with known pilot pattern:
% $M(l,k) = \sum_{l'=0}^{L-1}\sum_{k'=0}^{K-1} R(l+l',k+k') \cdot P^*(l',k')$
% Frame start indices detected from correlation peaks
Plk = get_drm_pilot_frame(iModeEst,stDRM.occupancy);
Mlk = xcorr2(Rlk,Plk);
Mlk(1:iNOfSymbolsPerFrame-1,:) = [];

if SwitchDemoSync
    figure(201);
    subplot(2,1,1)
    mesh([0:2*iNfft-2]-iNfft,0:size(Mlk,1)-1,abs(Mlk))

    xlabel('k Subchannel')
    ylabel('l Symbols')
    zlabel('|M(l,k)|')

    subplot(2,1,2)
    plot(0:size(Mlk,1)-1,abs(Mlk(:,iNfft)))
    xlabel('l Symbols')
    ylabel('|M(l,0)|')
end

iNOfFrames = floor(size(Mlk,1)/iNOfSymbolsPerFrame);

[~ , viFrameStart] = maxk(abs(Mlk(:,iNfft)),iNOfFrames);

viFrameStart = sort(viFrameStart);
viFrameStart(end)=[];

viFrameSymbols = 0:iNOfSymbolsPerFrame-1;

for iFrame = 1:length(viFrameStart)
    stRlk{iFrame} = Rlk(viFrameStart(iFrame)+viFrameSymbols,:);
end

% Run first Frame
Rlk = vertcat(stRlk{1:iNofFramesNeeded});
iNOfSymbolsTotal = iNofFramesNeeded * iNOfSymbolsPerFrame;

%% Fine Synchronization
% Corrects residual phase errors using pilot tones
% For each symbol l:
% 1. Calculate phase difference: $\phi(l,k) = \angle(R(l,k) \cdot P^*(l,k))$
% 2. Linear regression: $\phi(l,k) \approx ak + b$
% 3. Correction: $R_{corr}(l,k) = R(l,k) \cdot e^{-j(ak+b)}$
Plk = get_drm_pilot_frame(iModeEst,stDRM.occupancy);

% Apply fine synchronization per frame
for iFrame = 1:iNofFramesNeeded
    icurFrameStart = (iFrame-1)*iNOfSymbolsPerFrame+1;
    icurFrameEnd = iFrame*iNOfSymbolsPerFrame;
    Rlk(icurFrameStart:icurFrameEnd,:) = fine_sync(Rlk(icurFrameStart:icurFrameEnd,:), Plk, iNfft, SwitchDemoSync);
end

%% Channel Estimation and Equalization
% Channel estimation at pilot positions: $H(k_p) = \frac{R(k_p)}{P(k_p)}$
% Interpolation methods:
% - Spline: Cubic spline interpolation between pilots
% - Wiener: MMSE estimation based on channel statistics
%   $\hat{H}(k) = \mathbf{w}^H \mathbf{H_p}$ where $\mathbf{w} = \mathbf{R}^{-1}\mathbf{r}$
%   $\mathbf{R} = E[\mathbf{H_p}\mathbf{H_p}^H]$, $\mathbf{r} = E[\mathbf{H_p}H^*(k)]$
% Equalization: $S_{est}(k) = \frac{R(k)}{\hat{H}(k)}$
cInterpolater = 'Wiener';  % 'Spline' or 'Wiener'
for iFrame = 1:iNofFramesNeeded
    icurFrameStart = (iFrame-1)*iNOfSymbolsPerFrame+1;
    icurFrameEnd = iFrame*iNOfSymbolsPerFrame;
    Rlk(icurFrameStart:icurFrameEnd,:) = channel_estimation_equalization(Rlk(icurFrameStart:icurFrameEnd,:), Plk, stOFDM.iNfft, stOFDM.iNg, iModeEst, stDRM.occupancy, cInterpolater, SwitchDemoSync);
end

%% Real Bit Error Rate Calculation
% TODO: Implement the BER calculation
% BER is the ratio of bit errors to total bits: $BER = \frac{N_{errors}}{N_{total}}$
% [numberErrors, BER_real] = biterr(dataInBits, dataOutBits);

%% Reconstruct Image
% Demaps received symbols to bits, deinterleaves, and reconstructs the image
% Image data extraction follows reverse of transmission process
[reconstructed_image, Received_call_sign] = reconstruct_drm_image(Rlk, stDRM, M, image_size, call_sign);

% Display the call sign
fprintf('Received call sign: %s\n', char(Received_call_sign));

%% Graphical Output
% Shows comparison between transmitted and received data
SlkTemp = repmat(get_drm_data_template_frame(stDRM.mode, stDRM.occupancy), iNofFramesNeeded, 1);
figure(1)
subplot(2,3,1)
imagesc((1:stOFDM.iNfft)-get_drm_dc_position(stDRM.mode,stDRM.occupancy),0:size(Slk,1)-1,abs(Slk))
xlabel('Subchannel k')
ylabel('Symbol l')
title('Transmit Frame |Slk|')

subplot(2,3,2)
plot(Slk(SlkTemp ==1),'r.')
grid
fLimit = max(max(abs(Slk(SlkTemp ==1))));
axis square
axis([-fLimit fLimit -fLimit fLimit])
xlabel('I')
ylabel('Q')
title('Transmit Constellation')

subplot(2,3,3)
image = imread(image_path);
image = rgb2gray(image);

imshow(uint8(image));

subplot(2,3,4)
imagesc((1:stOFDM.iNfft)-get_drm_dc_position(stDRM.mode,stDRM.occupancy),0:size(Rlk,1)-1,abs(Rlk))
xlabel('Subchannel k')
ylabel('Symbol l')
title('Receive Frame |Rlk|')

subplot(2,3,5)
plot(Rlk(SlkTemp ==1),'r.')
grid
fLimit = max(max(abs(Rlk(SlkTemp ==1))));
axis square
axis([-fLimit fLimit -fLimit fLimit])
xlabel('I')
ylabel('Q')
title('Receive Constellation')

subplot(2,3,6)
imshow(uint8(reconstructed_image));


if SwitchDemoSync
    figure(2)
    subplot(3,2,[1 3 5])
    mfcPhase = angle(conj(Slk).*Rlk);
    imagesc((1:stOFDM.iNfft)-get_drm_dc_position(stDRM.mode,stDRM.occupancy),0:iNOfSymbolsTotal-1,mfcPhase)
    xlabel('Subchannel k')
    ylabel('Symbol l')
    title('Phase Difference')

    fpilot_position = get_drm_fpilot_position(stDRM.mode)+get_drm_dc_position(stDRM.mode,stDRM.occupancy);
    dc_position = get_drm_dc_position(stDRM.mode,stDRM.occupancy);
    subplot(3,2,2)
    plot(0:iNOfSymbolsTotal-1,mfcPhase(:,fpilot_position(1)))
    xlabel('Symbol l')
    ylabel('Phase Difference in rad')
    title(['1st fPilot / Subchannel:  ' num2str(fpilot_position(1)-dc_position)])
    grid


    subplot(3,2,4)
    plot(0:iNOfSymbolsTotal-1,mfcPhase(:,fpilot_position(2)))
    xlabel('Symbol l')
    ylabel('Phase Difference in rad')
    title(['2nd fPilot / Subchannel:  ' num2str(fpilot_position(2)-dc_position)])
    grid


    subplot(3,2,6)
    plot(0:iNOfSymbolsTotal-1,mfcPhase(:,fpilot_position(3)))
    xlabel('Symbol l')
    ylabel('Phase Difference in rad')
    title(['3rd fPilot / Subchannel:  ' num2str(fpilot_position(3)-dc_position)])
    grid
end
