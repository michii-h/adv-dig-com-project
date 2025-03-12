function [iModeEst, iNfft, iNg, iNs, iNOfSymbolsPerFrame] = detect_robustness_mode(vfcReceiveSignal, iOcc, SwitchDemoSync)
% DETECT_ROBUSTNESS_MODE Detects the DRM robustness mode from received signal
%
% Inputs:
%   vfcReceiveSignal - The received signal vector
%   iOcc - DRM occupancy parameter
%   SwitchDemoSync - Flag for enabling demo synchronization plots
%
% Outputs:
%   iModeEst - Estimated robustness mode (1=A, 2=B, 3=C, 4=D)
%   iNfft - FFT length for the detected mode
%   iNg - Guard interval length
%   iNs - Complete symbol length (iNfft + iNg)
%   iNOfSymbolsPerFrame - Number of symbols per DRM frame

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

% Determine OFDM Parameters for estimated Robustnes Mode
% FFT length
iNfft = get_drm_n_useful(iModeEst, iOcc);
% Guard Intervall Length
iNg = get_drm_n_guard(iModeEst, iOcc);
% Complete Symbol length
iNs = iNfft + iNg;
% Number of Symols per DRM Frame
iNOfSymbolsPerFrame = get_drm_symbols_per_frame(iModeEst);

strModes = ['A','B','C','D'];
fprintf('Robustness Mode: %s\n', strModes(iModeEst));
end
