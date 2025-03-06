%% OFDM Development Enviroment
start_up;

SwitchDemoSync = false; % additional plots for Demo

%% DRM Initialization
% stDRM Mode
stDRM.mode = 2; % Corresponds to Mode B
% Spectrum Occupancy
stDRM.occupancy = 3;

%% stOFDM Initialization
% FFT length
stOFDM.iNfft = get_drm_n_useful(stDRM.mode,stDRM.occupancy);
% Guard Intervall Length
stOFDM.iNg = get_drm_n_guard(stDRM.mode,stDRM.occupancy);
% Complete Symbol length
stOFDM.iNs = stOFDM.iNfft + stOFDM.iNg;

%% Generate DRM Frame
% Number of Symbols per frame
iNOfSymbols = get_drm_symbols_per_frame(stDRM.mode);

% Initialize Frame Matrix
Slk = get_drm_data_template_frame(stDRM.mode, stDRM.occupancy);

% Load Image
image = imread('image.png');
%image = image(1:20,1:20,:);
image_size = size(image);

viImage = image(:)';

% Reconstruct the image from the vector
% reconstructed_image = reshape(vector, size(image));

% Concat Call Sign
call_sign = uint8(double('DL0FHR'));
viData = [call_sign, viImage];

% mData = reshape( , 256);

% viData -> n x 256
%viData = reshape

%binaryData = de2bi(viData, 8, 'left-msb'); % convert to binary (8 bits per integer)
%binaryData = reshape(binaryData8.', [], 4); % reshape to 4 bits per row

% Datamapping: M-QAM
M = 256;     % 4 bits per row -> 16 possible values
%viDlk = qammod(binaryData,M,UnitAveragePower=true, InputType='bit');
viDlk = qammod(viData,M);
viDlk = reshape(viDlk.',1,[]); % concat rows one after another

% Set Data in Slk
vSlk = reshape(Slk',1,[]); % concat rows one after another

% Preallocate viDataPadded
DataPerFrame = sum(vSlk);

% Pad Data with Zeros at Pilot Positions
slkCtr = 1;
dataCtr = 1;
resultCtr = 1;
iNofFramesNeeded = ceil( size(viDlk, 2) / DataPerFrame );
viDataPadded = zeros(1, iNofFramesNeeded);
while dataCtr <= size(viDlk, 2)
    if vSlk(slkCtr) == 1
        viDataPadded(resultCtr) = viDlk(dataCtr);
        dataCtr = dataCtr + 1;
    else
        % already zero
    end
    slkCtr = mod(slkCtr, 15*256) + 1; % 15*256 = len of Slk
    resultCtr = resultCtr + 1;
end

% Padding at end for full DRM frames
nRows = ceil(size(viDataPadded, 2) / 256);
nRowsPad = 15 - mod(nRows, 15);
nRows = nRows + nRowsPad;

zData = zeros(1, nRows * 256);
zData(1:size(viDataPadded,2)) = viDataPadded;

viDataPadded = reshape(zData, 256, [])';
Dlk = viDataPadded;

% Slk(Slk == 1) = Dlk(Slk == 1);
Slk = Dlk;

% Generate Pilots
Plk = get_drm_pilot_frame(stDRM.mode,stDRM.occupancy);

% Duplicate Pilots to length of Slk
iNOfFrames = nRows / iNOfSymbols;
Plk = repmat(Plk,[iNOfFrames 1]);

% Set Pilots in Slk
Slk(Plk ~= 0) = Plk(Plk ~= 0);

%% OFDM Modulator
% IFFT
SlkTemp = ifft(fftshift(Slk,2),stOFDM.iNfft,2);

% Add Cyclic Prefix
SlkTemp = [SlkTemp(:,end-stOFDM.iNg+1:end) SlkTemp];

% Parallel to Serial Conversion
SlkTemp = SlkTemp.';
vfcTransmitSignal = SlkTemp(:);

% Repeat iG Frames
iG = 3;
vfcTransmitSignal = repmat(vfcTransmitSignal,iG,1);

%% Initialize Satellite Communication specific parameters for QO-100 - Adalm Pluto
% Link budget calculation etc
stSat = init_qo100_params();

%% Channel
% Switch Channels
iSwitchChannel = 2;

switch iSwitchChannel

    case 0 % ideal channel
        fprintf('Using ideal channel...\n');
        vfcReceiveSignal = vfcTransmitSignal;

    case 1 % simulated channel
        fprintf('Using simulated channel...\n');

        stChannel = initChannel();
        vfcReceiveSignal = channel_sim(vfcTransmitSignal,stChannel);

    case 2 % Simulate Satellite Communication for Q0-100
       vfcReceiveSignal = simulate_qo100_channel(vfcTransmitSignal, stSat);

    case 3 % Use Adalm Pluto
        fprintf('Using Adalm Pluto...\n');
        stAdalmPluto = initSDR(stSat.sampleRate, stSat.fc, stSat.adalm_txGain, stSat.adalm_rxGain);
        % Pass oversampling factor if needed
        if isfield(stSat, 'oversampling_factor')
            stAdalmPluto.oversampling_factor = stSat.oversampling_factor;
        end
        vfcReceiveSignal = LoopbackAdalmPluto(vfcTransmitSignal, stAdalmPluto);

end

%% Detect Robustness Mode
% force vfcReceiveSignal to be a column vector
vfcReceiveSignal = vfcReceiveSignal(:);

if SwitchDemoSync
    figure(100)
    plot(-length(vfcReceiveSignal)+1:length(vfcReceiveSignal)-1,abs(xcorr((vfcReceiveSignal))))
    xlabel('\Delta k')
    ylabel('r_{xx}(\Delta k)')
end

% Detect Robustness by autocorrelation at dk = [288 256 176 112]
% Autocorrelation matrix at dk corresponding to FFT lengths
mX(1,:) = circshift(vfcReceiveSignal,[288 0])';
mX(2,:) = circshift(vfcReceiveSignal,[256 0])';
mX(3,:) = circshift(vfcReceiveSignal,[176 0])';
mX(4,:) = circshift(vfcReceiveSignal,[112 0])';

% Calculate Robustnes Mode
[~, iModeEst] = max(abs(mX*vfcReceiveSignal));
iOcc = 3;
% Determine OFDM Parameters for estimated Robustnes Mode
% FFT length
iNfft = get_drm_n_useful(iModeEst,3); % Occupancy fixed at 3
% Guard Intervall Length
iNg = get_drm_n_guard(iModeEst,3);% Occupancy fixed at 3
% Complete Symbol length
iNs = iNfft + iNg;
% Number of Symols per DRM Frame
iNOfSymbolsPerFrame = get_drm_symbols_per_frame(iModeEst);

strModes = ['A','B','C','D'];
fprintf('Robustness Mode: %s\n', strModes(iModeEst));

%% Synchronization
% Number of Symbols in CaptureBuffer
iNOfSymbols = floor(length(vfcReceiveSignal)/iNs);

% Timing Recovery with Cyclic Prefix
vIndexGI = [0:iNg-1].';
R = zeros(1,iNs);

for l = 0:iNOfSymbols-2
    for k = 1:iNs
        R(k) = R(k) + ...
            sum(conj(vfcReceiveSignal(k+vIndexGI+l*iNs)).*vfcReceiveSignal(k+vIndexGI+iNfft+l*iNs));
    end
end

[~ ,iStartSample] = max(abs(R));
dOmegaEst = angle(R(iStartSample))/iNfft;

% Compensate Frequency Offset
vfcPhaser = exp(-1j*dOmegaEst*[0:length(vfcReceiveSignal)-1]);
vfcPhaser = vfcPhaser(:);

% Adjust to Start Sample
vfcReceiveSignal = vfcReceiveSignal(iStartSample:end);
iNOfSymbols = floor(length(vfcReceiveSignal)/iNs);
vfcReceiveSignal = vfcReceiveSignal(1:iNOfSymbols*iNs);

if SwitchDemoSync
    figure(101)
    plot([0:iNs-1],abs(R))
    xlabel('k')
    ylabel('R_{xy}(k)')
    hold on
    plot(iStartSample, abs(R(iStartSample)),'ro')
    hold off

    % figure;
    % plot(abs(sum(reshape(filter(ones(1,iNg),1,vfcReceiveSignal.*conj(circshift(vfcReceiveSignal,[iNfft 0]))),iNs,iNOfSymbols),2)))
end

%% OFDM Demodulator
% Serial to Parallel Conversion
RlkTemp = reshape(vfcReceiveSignal,stOFDM.iNs,[]).';

% Remove Cyclic Prefix
RlkTemp = RlkTemp(:,stOFDM.iNg+1:end);

% FFT
Rlk = fftshift(fft(RlkTemp,stOFDM.iNfft,2),2);

%% Frame Detection
% Generate Pilots
Plk = get_drm_pilot_frame(iModeEst,3);
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

% Run for full image
Rlk = vertcat(stRlk{1:end});
iNOfSymbols = iNofFramesNeeded * iNOfSymbolsPerFrame;
Rlk = Rlk(1:iNOfSymbols, :);

%% Fine Synchronization
Plk = get_drm_pilot_frame(iModeEst,iOcc);
% Plk = repmat(Plk,[iNofFramesNeeded 1]);

% Call the fine synchronization function for each frame to correct phase errors
for iFrame = 1:iNOfFramesNeeded
    icurFrameStart = (iFrame-1)*iNOfSymbolsPerFrame+1;
    icurFrameEnd = iFrame*iNOfSymbolsPerFrame;
    Rlk(icurFrameStart:icurFrameEnd,:) = fine_sync(Rlk(icurFrameStart:icurFrameEnd,:), Plk, iNfft, SwitchDemoSync);
end

%% Channel Estimation and Equalization
cInterpolater = 'Wiener';  % 'Spline' or 'Wiener'
Rlk = channel_estimation_equalization(Rlk, Plk, stOFDM.iNfft, stOFDM.iNg, iModeEst, iOcc, cInterpolater, SwitchDemoSync);

%% Reconstruct Image

SlkTemp = repmat(get_drm_data_template_frame(stDRM.mode, stDRM.occupancy), size(Rlk,1)/15, 1);

% Extrahiere Nutzdaten

% Initialisiere eine leere Liste für extrahierte Werte
viReceivedDlk = [];

% Gehe jede Zeile von Rlk durch und extrahiere die entsprechenden Werte
for iRow = 1:size(Rlk, 1)
    % Extrahiere nur die Werte in der aktuellen Zeile, wo SlkTemp == 1 ist
    viReceivedDlk = [viReceivedDlk, Rlk(iRow, SlkTemp(iRow, :) == 1)];
end

% QAM-Demodulation
binaryDataReceived = qamdemod(viReceivedDlk', M);

% Reshape zu einer Spaltenstruktur
%binaryDataReceived = reshape(binaryDataReceived, [], 8); % 8 Bits pro Integer

% Umwandlung in Integer-Daten
%viDataReceived = bi2de(binaryDataReceived, 8,'left-msb');

% Call-Sign entfernen (erste 6 Zeichen)
viImageReceived = binaryDataReceived(7:prod(image_size)+6);

% Zurück in 2D-Matrix mit urspr?nglicher Bildgr??e
reconstructed_image = reshape(viImageReceived, image_size);

% Anzeige des rekonstruierten Bildes
figure(140)
imshow(uint8(reconstructed_image));
%% Graphical Output
SlkTemp = repmat(get_drm_data_template_frame(stDRM.mode, stDRM.occupancy), iNofFramesNeeded, 1);
figure(1)
subplot(2,2,1)
imagesc((1:stOFDM.iNfft)-get_drm_dc_position(stDRM.mode,stDRM.occupancy),0:size(Slk,1)-1,abs(Slk))
xlabel('Subchannel k')
ylabel('Symbol l')
title('Transmit Frame |Slk|')

subplot(2,2,2)
plot(Slk(SlkTemp ==1),'r.')
grid
axis square
xlabel('I')
ylabel('Q')
title('Transmit Constellation')

subplot(2,2,3)
imagesc((1:stOFDM.iNfft)-get_drm_dc_position(stDRM.mode,stDRM.occupancy),0:size(Rlk,1)-1,abs(Rlk))
xlabel('Subchannel k')
ylabel('Symbol l')
title('Receive Frame |Rlk|')

subplot(2,2,4)
plot(Rlk(SlkTemp ==1),'r.')
grid
fLimit = max(max(abs(Rlk(SlkTemp ==1))));
axis square
axis([-fLimit fLimit -fLimit fLimit])

xlabel('I')
ylabel('Q')
title('Receive Constellation')


if SwitchDemoSync
    figure(2)
    subplot(3,2,[1 3 5])
    mfcPhase = angle(conj(Slk).*Rlk);
    imagesc((1:stOFDM.iNfft)-get_drm_dc_position(stDRM.mode,stDRM.occupancy),0:iNOfSymbols-1,mfcPhase)
    xlabel('Subchannel k')
    ylabel('Symbol l')
    title('Phase Difference')

    fpilot_position = get_drm_fpilot_position(stDRM.mode)+get_drm_dc_position(stDRM.mode,stDRM.occupancy);
    dc_position = get_drm_dc_position(stDRM.mode,stDRM.occupancy);
    subplot(3,2,2)
    plot(0:iNOfSymbols-1,mfcPhase(:,fpilot_position(1)))
    xlabel('Symbol l')
    ylabel('Phase Difference in rad')
    title(['1st fPilot / Subchannel:  ' num2str(fpilot_position(1)-dc_position)])
    grid


    subplot(3,2,4)
    plot(0:iNOfSymbols-1,mfcPhase(:,fpilot_position(2)))
    xlabel('Symbol l')
    ylabel('Phase Difference in rad')
    title(['2nd fPilot / Subchannel:  ' num2str(fpilot_position(2)-dc_position)])
    grid


    subplot(3,2,6)
    plot(0:iNOfSymbols-1,mfcPhase(:,fpilot_position(3)))
    xlabel('Symbol l')
    ylabel('Phase Difference in rad')
    title(['3rd fPilot / Subchannel:  ' num2str(fpilot_position(3)-dc_position)])
    grid
end
