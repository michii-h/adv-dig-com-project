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
    stSat.nb_uplink_start = 2400.050;  % MHz
    stSat.nb_uplink_end = 2400.300;    % MHz
    stSat.nb_downlink_start = 10489.550; % MHz
    stSat.nb_downlink_end = 10489.800;   % MHz

    % Center frequencies (Hz) - calculated from transponder parameters
    stSat.uplinkFreq = (stSat.nb_uplink_start + stSat.nb_uplink_end) * 1e6 / 2;
    stSat.downlinkFreq = (stSat.nb_downlink_start + stSat.nb_downlink_end) * 1e6 / 2;
    stSat.transponderBW = 2.7e3;       % Transponder bandwidth (Hz)

    % SDR parameters
    stSat.sampleRate = 1e6;            % 1 MSPS for Adalm Pluto
    stSat.oversampling_factor = 10;    % The actual sample rate will be multiplied by this in the SDR
    stSat.adalm_rxGain = 60;           % Initial RX gain (dB)
    % For Adalm Pluto, TX gain range is -89.75 to 0 dB (0 is maximum power)
    stSat.adalm_minTxGain = -89.75;    % Minimum TX gain for Adalm Pluto
    stSat.adalm_maxTxGain = 0;         % Maximum TX gain for Adalm Pluto (0 is maximum power)
    stSat.adalm_bpsk_symbolRate = 2400; % Common symbol rate for NB digital modes

    % Link budget parameters
    stSat.eirp = 39;                   % Satellite EIRP (dBW)
    stSat.gt = -13;                    % G/T ratio (assumed performance at satellite)
    stSat.antennaDiam = 0.9;           % Ground antenna diameter (meters)
    stSat.targetSNR = 10;              % Target SNR (dB)
    stSat.slantRange = 35786e3;        % Geostationary orbit distance (m)

    % Calculate link budget
    [stSat.txPower, stSat.linkMargin, stSat.expectedSNR] = calculateLinkBudget(stSat);

    % Adjust Adalm Pluto TX gain based on calculated power
    % For Adalm Pluto, 0dB = max power, and negative values reduce power
    powerInDBm = 10*log10(stSat.txPower) + 30;  % Convert W to dBm

    % Normalize to Adalm Pluto's scale (-89.75 to 0 dB)
    % Higher powerInDBm means we want more power (closer to 0)
    normalizedTxGain = -max(0, 30 - powerInDBm);  % Simple scaling approach
    stSat.adalm_txGain = max(stSat.adalm_minTxGain, min(stSat.adalm_maxTxGain, normalizedTxGain));

    % Print results
    fprintf('Required transmit power: %.2f dBW (%.2f W)\n', 10*log10(stSat.txPower), stSat.txPower);
    fprintf('Required transmit power: %.2f dBm\n', powerInDBm);
    fprintf('Link margin:  %.2f dB\n', stSat.linkMargin);
    fprintf('Expected SNR: %.2f dB\n', stSat.expectedSNR);
    fprintf('Adalm Pluto TX gain setting: %.2f dB\n', stSat.adalm_txGain);
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

    % Antenna gain calculation (with efficiency included)
    antennaGain = 10*log10(((pi*stSat.antennaDiam)/lambda)^2 * antennaEfficiency);

    % Free space path loss (using actual slant range)
    fsl = 20*log10(4*pi*stSat.slantRange/lambda);

    % Cable loss calculation
    totalCableLoss = cableLength * cableLoss;

    % System noise calculation
    noiseTemp = 290;        % Ambient temperature (K)
    systemTemp = noiseTemp / 10^(stSat.gt/10);
    noiseFloor = 10*log10(stSat.k * systemTemp * stSat.transponderBW);

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
    fprintf(' Required TX power: %.2f dBW (%.4f W)\n', requiredTxPowerDB, txPower);
end
