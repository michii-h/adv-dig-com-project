function [stSat] = main(antenna_diameter, tx_power_W)
% INIT_QO100_PARAMS Initialize parameters for QO-100 satellite
%   This function sets up all necessary parameters for QO-100 satellite link
%   using MATLAB's Satellite Communications Toolbox
%
%   Parameters:
%   antenna_diameter - Diameter of ground station antenna in meters
%   tx_power_W - Transmit power in Watts
%
%   Returns:
%   stSat - Structure with satellite parameters containing:
%     .name             - Satellite name
%     .sat_longitude    - Satellite longitude (degrees)
%     .orbit_altitude   - Altitude of orbit (meters)
%     .orbit_inclination- Inclination of orbit (degrees)
%     .uplink_freq      - Uplink frequency (Hz)
%     .downlink_freq    - Downlink frequency (Hz)
%     .frequency        - Reference frequency for Doppler (Hz)
%     .bandwidth        - Transponder bandwidth (Hz)
%     .transponder_power- Satellite transponder power (W)
%     .noise_figure     - Noise figure (dB)
%     .gs_latitude      - Ground station latitude (degrees)
%     .gs_longitude     - Ground station longitude (degrees)
%     .gs_altitude      - Ground station altitude (meters)
%     .range            - Range to satellite (km)
%     .az               - Azimuth to satellite (degrees)
%     .el               - Elevation to satellite (degrees)
%     .antenna_diameter - Ground station antenna diameter (m)
%     .antenna_efficiency - Antenna efficiency factor
%     .tx_gain_dB       - Transmit antenna gain (dBi)
%     .rx_gain_dB       - Receive antenna gain (dBi)
%     .tx_power         - Transmit power (W)
%     .tx_power_dBm     - Transmit power (dBm)
%     .tx_power_dBW     - Transmit power (dBW)
%     .path_loss_up     - Uplink path loss (dB)
%     .path_loss_down   - Downlink path loss (dB)
%     .uplink_budget    - Complete uplink budget with margins
%     .downlink_budget  - Complete downlink budget with margins
%     .SNR              - Signal-to-Noise ratio (dB)
%     .CNR              - Carrier-to-Noise ratio (dB)

% Input validation
validateattributes(antenna_diameter, {'numeric'}, {'scalar', 'positive', 'finite'}, 'init_qo100_params', 'antenna_diameter');
validateattributes(tx_power_W, {'numeric'}, {'scalar', 'positive', 'finite'}, 'init_qo100_params', 'tx_power_W');

% Physical constants
SPEED_OF_LIGHT = physconst('LightSpeed'); % m/s
BOLTZMANN_CONSTANT = physconst('Boltzmann'); % J/K

% Initialize the satellite structure with basic parameters
stSat = initBasicParams(antenna_diameter, tx_power_W);

% Calculate satellite position and link parameters
stSat = calculateSatellitePos(stSat);
stSat = calculateAntennaGains(stSat, SPEED_OF_LIGHT);
stSat = calculateLinkBudget(stSat, SPEED_OF_LIGHT, BOLTZMANN_CONSTANT);
end

function stSat = initBasicParams(antenna_diameter, tx_power_W)
    % Initialize the satellite structure with basic parameters
    stSat = struct();

    % Fixed QO-100 parameters
    stSat.name = 'QO-100';
    stSat.sat_longitude = 25.8;         % Geostationary at 25.8°E
    stSat.orbit_altitude = 35786000;    % GEO altitude in meters
    stSat.orbit_inclination = 0.05;     % Small inclination in degrees (not perfectly geostationary)
    stSat.uplink_freq = 2400e6 + 1.5e6; % 2.4 GHz with 1.5 MHz offset for narrowband
    stSat.downlink_freq = 10489.5e6;    % 10.489 GHz
    stSat.frequency = stSat.uplink_freq; % For Doppler calculation
    stSat.bandwidth = 2.7e3;            % 2.7 kHz for narrowband digital modes
    stSat.transponder_power = 10;       % 10 W for narrowband transponder
    stSat.noise_figure = 2.5;           % 2.5 dB noise figure

    % Ground station parameters - TH-Rosenheim R-Bau
    stSat.gs_latitude = 47.86603555991787;   % latitude in degrees
    stSat.gs_longitude = 12.108404806448695; % longitude in degrees
    stSat.gs_altitude = 447;                 % altitude in meters
    stSat.rain_rate = 25;                    % Average rain rate in mm/h (typical for Bavaria)

    % Antenna parameters
    stSat.antenna_diameter = antenna_diameter;
    stSat.antenna_efficiency = 0.5;     % 50% efficiency

    % Transmitter parameters
    stSat.tx_power = tx_power_W;
    stSat.tx_power_dBm = 10*log10(stSat.tx_power*1000);
    stSat.tx_power_dBW = 10*log10(stSat.tx_power);  % Convert W to dBW
    stSat.cable_loss = 3;               % 3 dB cable loss
    stSat.pointing_loss = 1;            % 1 dB pointing loss

    % System noise temperature
    stSat.system_temp_K = 290;          % System noise temperature in Kelvin

    % Satellite transponder parameters
    stSat.sat_rx_gain = 15;             % Satellite receive antenna gain in dBi
    stSat.sat_tx_gain = 25;             % Satellite transmit antenna gain in dBi
    stSat.sat_noise_figure = 2.0;       % Satellite receiver noise figure in dB
    stSat.sat_transponder_gain = 120;   % Satellite transponder gain in dB
end

function stSat = calculateSatellitePos(stSat)
    try
        % Create a satellite scenario
        startTime = datetime('now');
        stopTime = startTime + hours(1);
        sampleTime = 60;  % 60 seconds
        sc = satelliteScenario(startTime, stopTime, sampleTime);

        % Try to use TLE first, fall back to manual GEO setup
        try
            sat = satellite(sc, 'qo-100.tle');
        catch
            % If TLE file not found, create a GEO satellite manually
            sat = satellite(sc, 'SatelliteName', stSat.name, ...
                           'OrbitType', 'geostationary', ...
                           'Longitude', stSat.sat_longitude);
        end

        % Add the ground station
        gs = groundStation(sc, 'Name', 'TH-Rosenheim', ...
                       'Latitude', stSat.gs_latitude, ...
                       'Longitude', stSat.gs_longitude, ...
                       'Altitude', stSat.gs_altitude);

        % Calculate range, azimuth and elevation
        [azimuth, elevation, range] = aer(gs, sat);

        % Store the results in the structure
        stSat.range = range(1)/1000;  % Convert to km
        stSat.az = azimuth(1);
        stSat.el = elevation(1);
    catch e
        error('Error in satellite scenario calculation: %s\n', e.message);
    end
end

function stSat = calculateAntennaGains(stSat, c)
    % Calculate antenna gains using correct formulas
    % The gain of a parabolic antenna is: G = (η * π^2 * D^2) / λ^2
    % In dB: G_dB = 10*log10(η) + 20*log10(π*D/λ)

    % Uplink
    wavelength_up = c/stSat.uplink_freq;
    stSat.tx_gain_dB = 10*log10(stSat.antenna_efficiency) + 20*log10(pi*stSat.antenna_diameter/wavelength_up);

    % Downlink
    wavelength_down = c/stSat.downlink_freq;
    stSat.rx_gain_dB = 10*log10(stSat.antenna_efficiency) + 20*log10(pi*stSat.antenna_diameter/wavelength_down);
end

function stSat = calculateLinkBudget(stSat, c, k)
    % Calculate link budget using basic physics formulas

    % 1. Calculate Free Space Path Loss for both uplink and downlink
    % FSPL (dB) = 20*log10(d) + 20*log10(f) + 20*log10(4π/c)
    % where d is distance in meters, f is frequency in Hz

    d_m = stSat.range * 1000; % Convert km to m

    % Atmospheric attenuation estimates based on frequency and elevation
    % This is a simplified model; detailed models would use ITU recommendations
    atm_loss_up = calculateAtmosphericLoss(stSat.uplink_freq, stSat.el, stSat.rain_rate);
    atm_loss_down = calculateAtmosphericLoss(stSat.downlink_freq, stSat.el, stSat.rain_rate);

    % Calculate Free Space Path Loss (FSPL)
    stSat.path_loss_up = 20*log10(d_m) + 20*log10(stSat.uplink_freq) + 20*log10(4*pi/c) + atm_loss_up;
    stSat.path_loss_down = 20*log10(d_m) + 20*log10(stSat.downlink_freq) + 20*log10(4*pi/c) + atm_loss_down;

    % 2. Calculate Uplink Budget
    % EIRP = Tx Power (dBW) + Tx Antenna Gain (dBi) - Losses (dB)
    uplink_eirp = stSat.tx_power_dBW + stSat.tx_gain_dB - stSat.cable_loss - stSat.pointing_loss;

    % Power at satellite receiver = EIRP - path loss + satellite receive antenna gain
    uplink_received_power = uplink_eirp - stSat.path_loss_up + stSat.sat_rx_gain;

    % Noise power at satellite receiver
    Ts_sat = 290 * 10^(stSat.sat_noise_figure/10); % Satellite system noise temperature
    uplink_noise_power = 10*log10(k * Ts_sat * stSat.bandwidth);

    % Uplink CNR
    uplink_cnr = uplink_received_power - uplink_noise_power;

    % 3. Calculate Downlink Budget
    % Satellite EIRP (includes transponder gain applied to received signal)
    % We assume the transponder is operating in linear region
    sat_eirp = uplink_received_power + stSat.sat_transponder_gain + stSat.sat_tx_gain;

    % Power at ground receiver
    downlink_received_power = sat_eirp - stSat.path_loss_down + stSat.rx_gain_dB - stSat.cable_loss;

    % Ground station noise power
    Ts_gs = stSat.system_temp_K * 10^(stSat.noise_figure/10);
    downlink_noise_power = 10*log10(k * Ts_gs * stSat.bandwidth);

    % Downlink CNR
    downlink_cnr = downlink_received_power - downlink_noise_power;

    % 4. Overall link CNR (end-to-end)
    % For a bent-pipe transponder, the total CNR is:
    % 1/CNR_total = 1/CNR_uplink + 1/CNR_downlink
    uplink_cnr_linear = 10^(uplink_cnr/10);
    downlink_cnr_linear = 10^(downlink_cnr/10);
    total_cnr_linear = 1 / (1/uplink_cnr_linear + 1/downlink_cnr_linear);
    total_cnr_dB = 10*log10(total_cnr_linear);

    % Store the results
    stSat.uplink_budget = struct('EIRP', uplink_eirp, 'PathLoss', stSat.path_loss_up, ...
                               'ReceivedPower', uplink_received_power, ...
                               'NoisePower', uplink_noise_power, 'CNR', uplink_cnr);

    stSat.downlink_budget = struct('EIRP', sat_eirp, 'PathLoss', stSat.path_loss_down, ...
                                 'ReceivedPower', downlink_received_power, ...
                                 'NoisePower', downlink_noise_power, 'CNR', downlink_cnr);

    stSat.CNR = total_cnr_dB;
    stSat.SNR = total_cnr_dB; % For digital systems, CNR ≈ SNR with appropriate normalization
end

function atm_loss = calculateAtmosphericLoss(freq_Hz, elevation_deg, rain_rate_mm_hr)
    % Simplified atmospheric loss model based on frequency and elevation
    % This is a basic approximation; professional tools use ITU-R models

    freq_GHz = freq_Hz / 1e9;

    % Base atmospheric loss increases with frequency and decreases with elevation
    if freq_GHz < 5
        % Lower frequencies experience less atmospheric loss
        base_loss = 0.1 * freq_GHz;
    elseif freq_GHz < 10
        base_loss = 0.2 * freq_GHz;
    else
        % Higher frequencies experience more atmospheric loss
        base_loss = 0.3 * freq_GHz;
    end

    % Elevation correction - lower elevations have longer path through atmosphere
    elevation_factor = 1 / sin(deg2rad(max(elevation_deg, 5))); % Limit to minimum 5° elevation

    % Rain attenuation - simplified model
    % Increases with frequency and rain rate
    if rain_rate_mm_hr > 0
        if freq_GHz < 5
            rain_loss = 0.01 * rain_rate_mm_hr * freq_GHz;
        elseif freq_GHz < 10
            rain_loss = 0.05 * rain_rate_mm_hr * freq_GHz/5;
        else
            rain_loss = 0.1 * rain_rate_mm_hr * freq_GHz/10;
        end
    else
        rain_loss = 0;
    end

    % Total atmospheric loss
    atm_loss = base_loss * elevation_factor + rain_loss;
end
