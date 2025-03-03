function [stSat] = init_qo100_params(antenna_diameter, tx_power_W)
% INIT_QO100_PARAMS Initialize parameters for QO-100 satellite
%   This function sets up all necessary parameters for QO-100 satellite link
%   using MATLAB's Satellite Communications Toolbox
%
%   Parameters:
%   antenna_diameter - Diameter of ground station antenna in meters
%   tx_power_W - Transmit power in Watts
%
%   Returns:
%   stSat - Structure with satellite parameters

% Fixed QO-100 parameters
stSat.name = 'QO-100';
stSat.position = [25.8, 0, 35786];  % Geostationary at 25.8°E
stSat.uplink_freq = 2400e6 + 1.5e6; % 2.4 GHz with 1.5 MHz offset for narrowband
stSat.downlink_freq = 10489.5e6;    % 10.489 GHz
stSat.freq = stSat.uplink_freq;     % For Doppler calculation
stSat.bandwidth = 2.7e3;            % 2.7 kHz for narrowband digital modes
stSat.transponder_power = 10;       % 10 W for narrowband transponder
stSat.noise_figure = 2.5;           % 2.5 dB noise figure

% Ground station parameters - TH-Rosenheim R-Bau
stSat.gs_lat = 47.86603555991787;   % latitude
stSat.gs_lon = 12.108404806448695;  % longitude
stSat.gs_alt = 0.447;               % altitude in km

% Calculate satellite position relative to ground station using Satellite Communications Toolbox
% Create a satellite scenario
startTime = datetime('now');
stopTime = startTime + hours(1);
sampleTime = 60;  % 60 seconds
sc = satelliteScenario(startTime, stopTime, sampleTime);

% Add the QO-100 satellite as a geostationary satellite
tleFile = "q0-100.tle";
sat = satellite(sc, tleFile, 'Name', stSat.name);

% Add the ground station
gs = groundStation(sc, 'Name', 'TH-Rosenheim', 'Latitude', stSat.gs_lat, ...
                   'Longitude', stSat.gs_lon, 'Altitude', stSat.gs_alt*1000);  % Convert km to meters

% Get access data between satellite and ground station
ac = access(gs, sat);
accessIntervals = accessStatus(ac);

% Calculate range, azimuth and elevation
[range, azimuth, elevation] = rangeangle(gs, sat);

% Store the results in the structure
stSat.range = range(1)/1000;  % Convert to km
stSat.az = azimuth(1);
stSat.el = elevation(1);

% Calculate satellite velocity relative to Earth (zero for geostationary)
stSat.vel = [0, 0, 0];

% Antenna parameters
stSat.antenna_diameter = antenna_diameter;
stSat.antenna_efficiency = 0.5;     % 50% efficiency
stSat.antenna_gain = 20*log10((pi*stSat.antenna_diameter*stSat.uplink_freq/3e8)); % Gain in dBi

% Transmitter parameters
stSat.tx_power = tx_power_W;
stSat.tx_power_dBm = 10*log10(stSat.tx_power*1000);
stSat.cable_loss = 3;               % 3 dB cable loss
stSat.pointing_loss = 1;            % 1 dB pointing loss

% Link budget calculation
stSat.path_loss = 20*log10(4*pi*stSat.range*1000*stSat.uplink_freq/3e8);
stSat.EIRP = stSat.tx_power_dBm + stSat.antenna_gain - stSat.cable_loss - stSat.pointing_loss;
stSat.G_T = -13;                    % G/T of QO-100 satellite in dB/K
stSat.k = -228.6;                   % Boltzmann's constant in dBW/K/Hz
stSat.SNR = stSat.EIRP - stSat.path_loss + stSat.G_T - stSat.k - 10*log10(stSat.bandwidth);

% RF impairments for channel simulation
stSat.phase_noise = 0.01;           % Phase noise standard deviation
stSat.multipath = false;            % Minimal multipath for satellite link

% Calculate Doppler shift (for geostationary satellite this is very small)
stSat.doppler_shift = 0;  % Very close to zero for geostationary satellite

% Print link budget details
fprintf('QO-100 Link Budget Calculation:\n');
fprintf('  Uplink Frequency    = %.2f MHz\n', stSat.uplink_freq/1e6);
fprintf('  Antenna Diameter    = %.2f m\n', stSat.antenna_diameter);
fprintf('  Antenna Gain        = %.2f dBi\n', stSat.antenna_gain);
fprintf('  TX Power            = %.2f dBm (%.2f W)\n', stSat.tx_power_dBm, stSat.tx_power);
fprintf('  EIRP                = %.2f dBm\n', stSat.EIRP);
fprintf('  Path Loss           = %.2f dB\n', stSat.path_loss);
fprintf('  G/T of Satellite    = %.2f dB/K\n', stSat.G_T);
fprintf('  Bandwidth           = %.2f kHz\n', stSat.bandwidth/1e3);
fprintf('  Expected SNR        = %.2f dB\n', stSat.SNR);
fprintf('  Range to Satellite  = %.2f km\n', stSat.range);
fprintf('  Elevation Angle     = %.2f degrees\n', stSat.el);
fprintf('  Azimuth Angle       = %.2f degrees\n', stSat.az);

end
