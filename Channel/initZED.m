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
function stZedboard = initZedboard()

% Init Adalm Pluto        
stZedboard.fs = 12e3;     % Symbol rate
stZedboard.fc = 2.415e9; % Carrier center Frequency
stZedboard.TxGain = 0; % Gain, specified as a scalar from -89.75 to 0 dB with a resolution of 0.25 dB
stZedboard.RxGain = 0;  % Radio receiver gain in dB, specified as a scalar from -4 to 71

stZedboard.IPAddress='192.168.3.2';
stZedboard.TxChannelMapping=[1 2]; % Channel Mapping
stZedboard.RxChannelMapping=[1 2]; % Channel Mapping