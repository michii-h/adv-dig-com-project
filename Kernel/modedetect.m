function [ iMode ] = modedetect( x_baseband )
%% Detection Robustness Mode
% autocorrelation for shifts: 288, 256 , 176, 112 (Nfft of robustness
% modes)
% force column vector
x_baseband = x_baseband(:);
m_baseband = [ circshift(x_baseband,[288 0])' ; ...
               circshift(x_baseband,[256 0])' ; ...
               circshift(x_baseband,[176 0])' ; ...
               circshift(x_baseband,[112 0])' ];
vfcMetric =  m_baseband * x_baseband;

[value iMode] = max(abs(vfcMetric));

end

