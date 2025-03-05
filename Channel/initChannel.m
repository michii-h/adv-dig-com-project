% *************************************************************************
% Copyright: (c) 2022 Technische Hochschule Rosenheim  All rights reserved.
%                Hochschulstr. 1
%                D-83024 Rosenheim
% MODULE:        Initialisierung Channel Struktur
% PROJECT:       The Single Carrier Project
% COMPILER:      Matlab R2022a
% AUTHOR:        Markus Stichler
% HISTORY:		 15.06.2022
% *************************************************************************
% Initialisierung Channel
% ************************************************************************
function stChannel = initChannel()

% Initialize Channel settings
stChannel.fPhaseOffset          = 0;      % Phase Offset [rad]
stChannel.fFreqOffset           = 0.003  ;    % Frequency Offset [fs/iNfft]
stChannel.fClockOffsetPpm       = 0;      % Clock Offset [ppm]
stChannel.fSampleOffset         = 0.01;      % Sample Offset [bin]
stChannel.fOmegaOffset          = 0;
stChannel.vfImpulseResponse     = [1];      % Impulse response
stChannel.fSNRdB                = 1000;    % SNR [dB]
stChannel.fIqOffsetdB           = -100;   % IQ Offset [dB]
stChannel.fGainImbalancedB      = 0;      % Gain Imbalance [dB]
stChannel.fQuadratureError      = 0;      % Quadratur Error [rad]
stChannel.fPhaseDevMag          = 0;      % Phase Deviation Magnitude [rad]
stChannel.fPhaseDevPeriodLength = Inf;    % Phase Deviation Period Length [Samples]