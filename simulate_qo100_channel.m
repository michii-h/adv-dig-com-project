function [vfcReceiveSignal] = simulate_qo100_channel(vfcTransmitSignal, stSat, stOFDM)
% SIMULATE_QO100_CHANNEL Simulates the channel effects of the QO-100 satellite
%   This function simulates the effects of the QO-100 satellite channel on the
%   transmitted signal using the Satellite Communications Toolbox
%
%   Parameters:
%   vfcTransmitSignal - Complex baseband transmitted signal
%   stSat - Structure with satellite parameters including:
%       .orbit_altitude - Satellite altitude in meters (35786000 for GEO)
%       .orbit_inclination - Inclination in degrees (typically ~0 for GEO)
%       .frequency - Carrier frequency in Hz
%       .gs_latitude - Ground station latitude in degrees
%       .gs_longitude - Ground station longitude in degrees
%       .gs_altitude - Ground station altitude in meters
%       .sat_longitude - Satellite longitude in degrees
%       .phase_noise - Phase noise standard deviation
%       .rain_rate - Ground station rain rate in mm/h
%       .multipath - Boolean flag for multipath simulation
%       .tx_power_dBW - Transmit power in dBW
%       .tx_gain_dB - Transmit antenna gain in dB
%       .rx_gain_dB - Receive antenna gain in dB
%       .system_temp_K - System noise temperature in Kelvin
%   stOFDM - Structure with OFDM parameters
%
%   Returns:
%   vfcReceiveSignal - Complex baseband received signal after channel

% Calculate slant range to satellite
slantRange = slantRangeCircularOrbit(stSat.orbit_altitude, ...
                                    stSat.gs_latitude, ...
                                    stSat.gs_longitude, ...
                                    stSat.gs_altitude, ...
                                    stSat.sat_longitude);

% Calculate Doppler shift (minimal for GEO satellites)
if isfield(stSat, 'orbit_inclination')
    % Only calculate Doppler if satellite isn't perfectly geostationary
    doppler = dopplerShiftCircularOrbit(stSat.orbit_altitude, ...
                                       stSat.orbit_inclination, ...
                                       stSat.frequency, ...
                                       stSat.gs_latitude, ...
                                       stSat.gs_longitude, ...
                                       stSat.gs_altitude, ...
                                       0); % Time = 0 for current position
else
    doppler = 0;
end

% Apply frequency scaling due to Doppler
if doppler ~= 0
    vfcPhaser = exp(1j * 2 * pi * doppler * [0:length(vfcTransmitSignal)-1]' / stOFDM.iNs);
    vfcTransmitSignal = vfcTransmitSignal .* vfcPhaser;
end

% Calculate propagation losses using ITU-R P.618 model
[totalLoss, ~, ~] = p618PropagationLosses(stSat.frequency, ...
                                        stSat.gs_latitude, ...
                                        stSat.rain_rate, ...
                                        slantRange/1000, ... % Convert to km
                                        stSat.gs_altitude/1000); % Convert to km

% Apply path loss and phase rotation
path_loss = db2mag(-totalLoss);
vfcTransmitSignal = vfcTransmitSignal * path_loss;

% Add transponder non-linearity (TWTA model)
% QO-100 has a non-linear transponder that causes AM/AM and AM/PM distortion
input_power = abs(vfcTransmitSignal).^2;
norm_power = input_power / max(input_power);

% AM/AM distortion (Saleh model)
alpha_a = 2.1587;
beta_a = 1.1517;
am_am = alpha_a * norm_power ./ (1 + beta_a * norm_power);

% AM/PM distortion
alpha_p = 4.0033;
beta_p = 9.1040;
am_pm = alpha_p * norm_power ./ (1 + beta_p * norm_power);

% Apply non-linear distortion
distortion = am_am .* exp(1j * am_pm);
distortion = distortion / max(abs(distortion)); % Normalize
vfcTransmitSignal = vfcTransmitSignal .* distortion;

% Add phase noise
phase_noise = stSat.phase_noise * randn(size(vfcTransmitSignal));
vfcTransmitSignal = vfcTransmitSignal .* exp(1j * phase_noise);

% Calculate CNR using the satellite link budget model
if isfield(stSat, 'tx_power_dBW') && isfield(stSat, 'tx_gain_dB') && ...
   isfield(stSat, 'rx_gain_dB') && isfield(stSat, 'system_temp_K')

    cnr = satelliteCNR('TransmitPower', stSat.tx_power_dBW, ...
                      'TransmitterGain', stSat.tx_gain_dB, ...
                      'ReceiverGain', stSat.rx_gain_dB, ...
                      'FrequencyMHz', stSat.frequency/1e6, ...
                      'SystemNoiseTemperature', stSat.system_temp_K, ...
                      'PropagationLoss', totalLoss);
    SNR = cnr;
else
    % Fallback to the provided SNR if parameters aren't available
    SNR = stSat.SNR;
end

% Add AWGN based on calculated SNR
vfcReceiveSignal = awgn(vfcTransmitSignal, SNR, 'measured');

% Add multipath effects (if any)
if isfield(stSat, 'multipath') && stSat.multipath
    % QO-100 has minimal multipath, but we can model some reflections
    reflections = [1 0.1 0.05];
    delays = [0 5 10]; % in samples

    % Create multipath channel
    h_mp = zeros(1, max(delays) + 1);
    for i = 1:length(reflections)
        h_mp(delays(i) + 1) = reflections(i);
    end

    % Apply multipath
    vfcReceiveSignal = conv(vfcReceiveSignal, h_mp, 'same');
end

end
