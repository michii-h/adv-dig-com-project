%% stOFDM Workshop: Lesson One
clc
close all
clear all
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

%% Generate Frame Template
% Number of Symbols per frame
iNOfSymbols = get_drm_symbols_per_frame(stDRM.mode);

% Initialize Frame Matrix
Slk = get_drm_data_template_frame(stDRM.mode, stDRM.occupancy);

% Datamapping: M-QAM
M = 4;
viRand = randi([0 M-1],iNOfSymbols,stOFDM.iNfft);
if M == 4
    Dlk = qammod(viRand,M,'UnitAveragePower',true);
else
    Dlk = qammod(viRand,M);
end

% Set Data in Slk
Slk(Slk == 1) = Dlk(Slk == 1);

% Generate Pilots
Plk = get_drm_pilot_frame(stDRM.mode,stDRM.occupancy);

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

%% Channel
% Switch Cahnnel
iSwitchChannel = 1;

switch iSwitchChannel

    case 0 % ideal channel
        vfcReceiveSignal = vfcTransmitSignal;

    case 1 % simulated channel
        stChannel = initChannel();
        vfcReceiveSignal = channel_sim(vfcTransmitSignal,stChannel,stOFDM);

    case 2 % capture with TenTec
        vfcReceiveSignal = captureDRM(2) ;

    case 3 % use Adalm Pluto for TX and RX
       stAdalmPluto = initSDR();
       vfcReceiveSignal = LoopbackAdalmPluto(vfcTransmitSignal,stAdalmPluto);

    case 4 % use ZedBoard for TX and RX
       stZedboard = initZED();
       vfcReceiveSignal = LoopbackZedBoard(vfcTransmitSignal,stZedboard);
       vfcReceiveSignal =vfcReceiveSignal(:,1);

end
%% Detect Robustness Mode


%% Synchronization


%% OFDM Demodulator
% Serial to Parallel Conversion
RlkTemp = reshape(vfcReceiveSignal,stOFDM.iNs,iNOfSymbols).';

% Remove Cyclic Prefix
RlkTemp = RlkTemp(:,stOFDM.iNg+1:end);

% FFT
Rlk = fftshift(fft(RlkTemp,stOFDM.iNfft,2),2);

%% Frame Detection


%% Fine Synchronization


%% Kanalschätzung und Entzerrung



%% Graphical Output

figure(1)
subplot(2,2,1)
imagesc((1:stOFDM.iNfft)-get_drm_dc_position(stDRM.mode,stDRM.occupancy),0:size(Slk,1)-1,abs(Slk))
xlabel('Subchannel k')
ylabel('Symbol l')
title('Transmit Frame |Slk|')

subplot(2,2,2)
plot(Slk(Plk ==0),'r.')
grid
axis square
xlabel('I')
ylabel('Q')
title('Transmit Constellation')

subplot(2,2,3)
imagesc((1:stOFDM.iNfft)-get_drm_dc_position(stDRM.mode,stDRM.occupancy),0:iNOfSymbols-1,abs(Rlk))
xlabel('Subchannel k')
ylabel('Symbol l')
title('Receive Frame |Rlk|')

subplot(2,2,4)
plot(Rlk(Plk ==0),'r.')
grid
fLimit = max(max(abs(Rlk(Plk ==0))));
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