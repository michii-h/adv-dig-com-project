%% QO-100 Parameters Initialization
% Initialize parameters for the QO-100 satellite communications
%
% Returns:
%   stSat - Structure containing all satellite and link parameters
%
% Example:
%   params = init_qo100_params();
function stSat = init_qo100_params()
    % Physical constants
    stSat.c = 3e8;                % Speed of light (m/s)
    stSat.k = 1.38e-23;           % Boltzmann constant (J/K)

    % Basic satellite parameters
    stSat.name = 'Qatar-OSCAR 100';

    % QO-100 Narrowband Transponder parameters (MHz)
    stSat.nb_uplink_start = 2400.150;  % MHz
    stSat.nb_uplink_end = 2400.245;    % MHz
    stSat.nb_downlink_start = 10489.650; % MHz
    stSat.nb_downlink_end = 10489.745;   % MHz

    % Center frequencies (Hz) - calculated from transponder parameters
    stSat.uplinkFreq = (stSat.nb_uplink_start + stSat.nb_uplink_end) * 1e6 / 2;
    stSat.downlinkFreq = (stSat.nb_downlink_start + stSat.nb_downlink_end) * 1e6 / 2;

    transponderShift = 8089.5e6;
    Lo = 9750e6;
    stSat.downlinkFreq = stSat.uplinkFreq + transponderShift - Lo;

    % SDR parameters
    stSat.fs = 2.7e3;                 % Sample rate for SDR (Hz)
    stSat.fc = stSat.uplinkFreq;      % Center frequency (Hz)
    % stSat.fcTx = stSat.uplinkFreq;    % TX center frequency (Hz)
    % stSat.fcRx = stSat.downlinkFreq;  % RX center frequency (Hz)
    stSat.oversampling_factor = 25;  % Oversampling factor
    stSat.bandwidth = 2.7e3;          % Bandwidth limit (Hz)

    % Adalm Pluto specific parameters
    stSat.rxGain = 50;                % RX gain (dB)
    stSat.minTxGain = -89.75;         % Minimum TX gain for Adalm Pluto
    stSat.maxTxGain = 0;              % Maximum TX gain (0 is maximum power)
    stSat.max_power_dbm = 7 + 50;     % Maximum power in dBm

    % Link budget parameters
    stSat.eirp = 39;                  % Satellite EIRP (dBW)
    stSat.gt = -13;                   % G/T ratio (assumed performance at satellite)
    stSat.antennaDiam = 0.9;          % Ground antenna diameter (meters)
    stSat.targetSNR = 25;             % Target SNR (dB)
    stSat.slantRange = 35786e3;       % Geostationary orbit distance (m)

    % Channel simulation parameters
    stSat.phaseNoiseVariance = 0.01;  % Phase noise variance
    stSat.freqOffset = 50;            % Frequency offset in Hz
    stSat.filterOrder = 64;           % Order of the transponder filter

    % Saleh model parameters for nonlinear distortion
    stSat.saleh.alpha_a = 2.1587;     % AM/AM parameter
    stSat.saleh.beta_a = 1.1517;      % AM/AM parameter
    stSat.saleh.alpha_p = 4.0033;     % AM/PM parameter
    stSat.saleh.beta_p = 9.1040;      % AM/PM parameter

    % Control flags
    stSat.skipNonlinear = true;       % Set to true to bypass nonlinear modeling
    stSat.skipFiltering = true;       % Set to true to bypass transponder filtering

    % Calculate link budget
    [stSat.txPower, stSat.linkMargin, stSat.expectedSNR] = calculateLinkBudget(stSat);

    % Convert required power to Adalm Pluto TX gain setting
    % For Adalm Pluto, 0dB = max_dBm, and negative values reduce power
    powerInDBm = 10*log10(stSat.txPower) + 30;  % Convert W to dBm

    % Assume Adalm Pluto max power is approximately 10 dBm at 0 dB gain
    % So to get our desired power, we need to set gain to (desired dBm - max_dBm)
    % But since gain can only be negative, we cap it at 0
    stSat.txGain = min(0, powerInDBm - stSat.max_power_dbm);

    % Ensure within valid range
    stSat.txGain = max(stSat.minTxGain, min(stSat.maxTxGain, stSat.txGain));

    % Store the actual output power we expect from the SDR with this gain
    stSat.expected_output_power_dbm = stSat.max_power_dbm + stSat.txGain;
    stSat.expected_output_power_w = 10^((stSat.expected_output_power_dbm - 30)/10);

    % Print results
    fprintf('Required transmit power: %.2f dBW (%.2f W)\n', 10*log10(stSat.txPower), stSat.txPower);
    fprintf('Required transmit power: %.2f dBm\n', powerInDBm);
    fprintf('Adalm Pluto TX gain setting: %.2f dB\n', stSat.txGain);
    fprintf('Expected actual output power: %.2f dBm (%.4f W)\n',...
            stSat.expected_output_power_dbm, stSat.expected_output_power_w);
    fprintf('Link margin:  %.2f dB\n', stSat.linkMargin);
    fprintf('Expected SNR: %.2f dB\n', stSat.expectedSNR);
    fprintf('Note: Adalm Pluto TX gain range is -89.75 to 0 dB (0 is maximum power)\n');

    % Add warning if link margin is negative
    if stSat.linkMargin < 0
        warning('Link margin is negative (%.2f dB). Communication may be unreliable.', stSat.linkMargin);
    end
end

function [txPower, margin, snr] = calculateLinkBudget(stSat)
    % Common calculations
    lambda = stSat.c / stSat.uplinkFreq;

    % Define parameters for both calculation methods
    antennaEfficiency = 0.5;      % 50% efficiency as used in example
    cableLength = 12;             % meters of cable
    cableLoss = 0.3;              % dB/m for Ecoflex 7mm

    % Antenna gain calculation
    antennaGain = 10*log10(((pi*stSat.antennaDiam)/lambda)^2 * antennaEfficiency);

    % Free space path loss
    fsl = 20*log10(4*pi*stSat.slantRange/lambda);

    % Cable loss calculation
    totalCableLoss = cableLength * cableLoss;

    % System noise calculation
    noiseTemp = 3;
    systemTemp = noiseTemp / 10^(stSat.gt/10);
    noiseFloor = 10*log10(stSat.k * systemTemp * stSat.bandwidth);

    % noiseFloor = 10*log10(stSat.k * stSat.transponderBW) + stSat.gt;

    % Calculate required received power at satellite
    requiredRxPower = noiseFloor + stSat.targetSNR;

    % Calculate required transmit power (dBW)
    requiredTxPowerDB = requiredRxPower + fsl + totalCableLoss - antennaGain;

    % Convert to linear power (Watts)
    txPower = 10^(requiredTxPowerDB/10);

    % Calculate actual SNR with this power
    rxPowerAtSat = 10*log10(txPower) + antennaGain - fsl - totalCableLoss;
    snr = rxPowerAtSat - noiseFloor;
    margin = snr - stSat.targetSNR;

    % ---- Print results ----
    fprintf('Link budget calculations:\n');
    fprintf(' Free space path loss: %.2f dB\n', fsl);
    fprintf(' Antenna gain: %.2f dBi\n', antennaGain);
    fprintf(' Cable loss: %.2f dB\n', totalCableLoss);
    fprintf(' Noise floor: %.2f dBW\n', noiseFloor);
end
