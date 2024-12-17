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
function stAdalmPluto = initSDR()

% Init Adalm Pluto        
stAdalmPluto.fs = 12e3;     % Symbol rate
stAdalmPluto.fc = 2.415e9; % Carrier center Frequency
stAdalmPluto.TxGain = -10; % Gain, specified as a scalar from -89.75 to 0 dB with a resolution of 0.25 dB
stAdalmPluto.RxGain = 50;  % Radio receiver gain in dB, specified as a scalar from -4 to 71