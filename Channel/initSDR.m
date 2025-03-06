% *************************************************************************
% Copyright: (c) 2022 Technische Hochschule Rosenheim  All rights reserved.
%                Hochschulstr. 1
%                D-83024 Rosenheim
% MODULE:        Initialisierung Software Defined Radio
% PROJECT:       The Single Carrier Project
% COMPILER:      Matlab R2022a
% AUTHOR:        Markus Stichler
% HISTORY:		 15.06.2022
% *************************************************************************
% Initialisierung SDR
% ************************************************************************
function stAdalmPluto = initSDR(stSat)

% Init Adalm Pluto
% stAdalmPluto.fs = 12e3;     % Symbol rate
% stAdalmPluto.fc = 2.415e9; % Carrier center Frequency
% stAdalmPluto.TxGain = -10; % Gain, specified as a scalar from -89.75 to 0 dB with a resolution of 0.25 dB
% stAdalmPluto.RxGain = 50;  % Radio receiver gain in dB, specified as a scalar from -4 to 71

stAdalmPluto.fs = stSat.sampleRate;     % Symbol rate
stAdalmPluto.fc = stSat.uplinkFreq;             % Carrier center Frequency
stAdalmPluto.TxGain = stSat.adalm_txGain; % Gain, specified as a scalar from -89.75 to 0 dB with a resolution of 0.25 dB
stAdalmPluto.RxGain = stSat.adalm_rxGain;  % Radio receiver gain in dB, specified as a scalar from -4 to 71

stAdalmPluto.max_power_dbm  = stSat.adalm_max_power_dbm;

transponderShift = 8089.5e6;
Lo = 9750e6;
stAdalmPluto.fcRx = stSat.uplinkFreq + transponderShift - Lo;

end
