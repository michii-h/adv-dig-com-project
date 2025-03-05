%% QO-100 Parameters Initialization
function stSat = init_qo100_params()
    % Basic satellite parameters
    stSat.name = 'Qatar-OSCAR 100';

    % QO-100 Narrowband Transponder parameters (MHz)
    stSat.nb_uplink_start = 2400.050;  % MHz
    stSat.nb_uplink_end = 2400.300;    % MHz
    stSat.nb_downlink_start = 10489.550; % MHz
    stSat.nb_downlink_end = 10489.800;   % MHz

    % Center frequencies (Hz)
    stSat.uplinkFreq = 2.4e9;          % Exact center freq
    stSat.downlinkFreq = 10.489675e9;  % Exact center freq
    stSat.transponderBW = 250e3;       % 250 kHz narrowband transponder

    % SDR parameters
    stSat.sampleRate = 1e6;            % 1 MSPS for Adalm Pluto
    stSat.adalm_txGain = 0;            % Initial TX gain (dB), adjust based on link budget
    stSat.adalm_rxGain = 60;           % Initial RX gain (dB)
    stSat.adalm_maxTxGain = 89;        % Maximum TX gain for Adalm Pluto
    stSat.adalm_bpsk_symbolRate = 2400; % Common symbol rate for NB digital modes

    % Link budget parameters
    stSat.eirp = 39;                   % Satellite EIRP (dBW)
    stSat.gt = 10;                     % G/T ratio (dB/K)
    stSat.antennaDiam = 0.9;           % Ground antenna diameter (meters)
    stSat.rainMargin = 3;              % Rain fade margin (dB)
    stSat.targetSNR = 10;              % Target SNR (dB)
    stSat.slantRange = 35786e3;        % Geostationary orbit distance (m)

    % Calculate link budget
    [stSat.txPower, stSat.linkMargin] = calculateLinkBudget(stSat);

    % Adjust Adalm Pluto TX gain based on calculated power
    stSat.adalm_txGain = - min(stSat.adalm_maxTxGain, 10*log10(stSat.txPower) + 30);

    fprintf('Required transmit power: %.2f dBW (%.2f W)\n', 10*log10(stSat.txPower), stSat.txPower);
    fprintf('Link margin: %.2f dB\n', stSat.linkMargin);
    fprintf('Recommended Adalm Pluto TX gain: %.2f dB\n', stSat.adalm_txGain);
end

function [txPower, margin] = calculateLinkBudget(stSat)
    % Constants
    c = 3e8;                % Speed of light (m/s)
    k = 1.38e-23;           % Boltzmann constant

    % Wavelength calculation
    lambda = c / stSat.uplinkFreq;

    % Free space path loss calculation
    fsl = 20*log10(4*pi*stSat.slantRange/lambda);

    % Antenna gain calculation
    efficiency = 0.55;      % Typical dish efficiency
    antennaGain = 10*log10(efficiency * (pi * stSat.antennaDiam / lambda)^2);

    % System noise calculation
    noiseTemp = 290;        % Ambient temperature (K)
    systemTemp = noiseTemp / 10^(stSat.gt/10);
    noiseFloor = 10*log10(k * systemTemp * stSat.sampleRate);

    % Link budget calculation
    rxPowerAtSat = stSat.eirp - fsl + antennaGain - stSat.rainMargin;
    margin = rxPowerAtSat - noiseFloor - stSat.targetSNR;

    % Calculate required transmit power in Watts
    txPower = 10^((stSat.eirp - antennaGain) / 10);

    % Log detailed calculations for debugging
    fprintf('Free space path loss: %.2f dB\n', fsl);
    fprintf('Antenna gain: %.2f dBi\n', antennaGain);
    fprintf('Noise floor: %.2f dBW\n', noiseFloor);
end
