%% OFDM Development Enviroment
start_up;

SwitchDemoSync = false; % additional plots for Demo

%% Transmission Parameters Initialization
% Create SatelliteLink object with appropriate parameters
tx_power = 7 + 30; % 7 dBm output power of Adalm Pluto + 30 dBm amplifier
link = SatelliteLink('tx_center_frequency', 2400.113 * 1e6, 'tx_gain', 0, 'rx_gain', 30, 'tx_power', tx_power);
link.useSameFrequency = true;
link.displayLinkBudget();

%% Bit Error Rate Calculation from expected SNR
ModOrder = 16;
BER_sim = berawgn(link.SNR_linear, 'qam', ModOrder);

fprintf('  Simulated BER: %.2e\n', BER_sim);

%% DRM Initialization
stDRM.mode = 4;      % Corresponds to Mode D
stDRM.occupancy = 3;

check_drm_bandwidth(link.baseStation.baseband_sample_rate, stDRM.mode, stDRM.occupancy);

%% stOFDM Initialization
% FFT length
stOFDM.iNfft = get_drm_n_useful(stDRM.mode,stDRM.occupancy);
% Guard Intervall Length
stOFDM.iNg = get_drm_n_guard(stDRM.mode,stDRM.occupancy);
% Complete Symbol length
stOFDM.iNs = stOFDM.iNfft + stOFDM.iNg;

%% Generate DRM Frame
% Call the function to generate the DRM frame
image_path = 'th-rosenheim-logo-colored-square.png';
call_sign = 'DL0FHR';
[Slk, M, image_size, iNofFramesNeeded, iNOfFrames] = generate_drm_frames(stDRM, stOFDM, image_path, call_sign);

%% OFDM Modulator
% IFFT
SlkTemp = ifft(fftshift(Slk,2),stOFDM.iNfft,2);

% Add Cyclic Prefix
SlkTemp = [SlkTemp(:,end-stOFDM.iNg+1:end) SlkTemp];

% Parallel to Serial Conversion
SlkTemp = SlkTemp.';
vfcTransmitSignal = SlkTemp(:);

%% Channel
% Switch Channels
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

        % Compensate fixed freq offset
        % Plk = get_drm_pilot_frame(stDRM.mode, stDRM.occupancy);
        % spectrumPlk = abs(fft(Plk(:)));
        % spectrumRx = abs(fft(vfcReceiveSignal));

        % comb = zeros(length(spectrumRx),1);
        % comb(round(750*(12000/link.baseStation.baseband_sample_rate))) = 10;
        % comb(round(2250*(12000/link.baseStation.baseband_sample_rate))) = 10;
        % comb(round(3000*(12000/link.baseStation.baseband_sample_rate))) = 10;

        % [~, offset] = max(xcorr(comb, spectrumRx));

        % freqShift = link.baseStation.baseband_sample_rate * offset / length(spectrumRx);

        freqShift = 282; % TODO: Edit this with real offset measured in plot

        delta_Omega = 2*pi*(freqShift/link.baseStation.baseband_sample_rate);
        vPhaser = exp(-1i * delta_Omega * [0:length(vfcReceiveSignal)-1])';

        vfcReceiveSignal = vfcReceiveSignal .* vPhaser;

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
[iModeEst, iNfft, iNg, iNs, iNOfSymbolsPerFrame] = detect_robustness_mode(vfcReceiveSignal, stDRM.occupancy, SwitchDemoSync);

%% Synchronization
vfcReceiveSignal = sync(vfcReceiveSignal, iNs, iNg, iNfft, SwitchDemoSync);

%% OFDM Demodulator
% Serial to Parallel Conversion
RlkTemp = reshape(vfcReceiveSignal,stOFDM.iNs,[]).';

% Remove Cyclic Prefix
RlkTemp = RlkTemp(:,stOFDM.iNg+1:end);

% FFT
Rlk = fftshift(fft(RlkTemp,stOFDM.iNfft,2),2);

%% Frame Detection
% Generate Pilots
Plk = get_drm_pilot_frame(iModeEst,stDRM.occupancy);
% Calculate Correlation Metrik
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
Plk = get_drm_pilot_frame(iModeEst,stDRM.occupancy);
% Plk = repmat(Plk,[iNofFramesNeeded 1]);

% Call the fine synchronization function for each frame to correct phase errors
for iFrame = 1:iNofFramesNeeded
    icurFrameStart = (iFrame-1)*iNOfSymbolsPerFrame+1;
    icurFrameEnd = iFrame*iNOfSymbolsPerFrame;
    Rlk(icurFrameStart:icurFrameEnd,:) = fine_sync(Rlk(icurFrameStart:icurFrameEnd,:), Plk, iNfft, SwitchDemoSync);
end

%% Channel Estimation and Equalization
cInterpolater = 'Wiener';  % 'Spline' or 'Wiener'
for iFrame = 1:iNofFramesNeeded
    icurFrameStart = (iFrame-1)*iNOfSymbolsPerFrame+1;
    icurFrameEnd = iFrame*iNOfSymbolsPerFrame;
    Rlk(icurFrameStart:icurFrameEnd,:) = channel_estimation_equalization(Rlk(icurFrameStart:icurFrameEnd,:), Plk, stOFDM.iNfft, stOFDM.iNg, iModeEst, stDRM.occupancy, cInterpolater, SwitchDemoSync);
end

%% Real Bit Error Rate Calculation
% TODO: Implement the BER calculation

% [numberErrors, BER_real] = biterr(dataInBits, dataOutBits);

%% Reconstruct Image
% Call the external function to reconstruct the image
[reconstructed_image, Received_call_sign] = reconstruct_drm_image(Rlk, stDRM, M, image_size, call_sign);

% Display the call sign
fprintf('Received call sign: %s\n', char(Received_call_sign));

%% Graphical Output
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
