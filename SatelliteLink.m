% 'tx_center_frequency', 2400.172 * 1e6, ...   % Uplink center frequency [Hz]
% 'baseband_sample_rate', 2.7e3, ...   % Baseband sample rate [Hz]
% 'tx_gain', 0, ...                    % TX Gain [dB]
% 'rx_gain', 50, ...                   % RX Gain [dB]

classdef SatelliteLink
    % This class calculates the link budget for a satellite link considering
    % base station and satellite parameters. The base station parameters
    % (TX/RX center frequency, baseband sample rate, TX/RX gain, output power,
    % cable losses, and antenna data) are grouped in a dedicated structure field 'baseStation'.

    properties
        % Base station parameters as a structure
        baseStation = struct( ...
            'tx_center_frequency', 2400.172 * 1e6, ...   % Uplink center frequency [Hz]
            'baseband_sample_rate', 2.7e3, ...  % Baseband sample rate [Hz]
            'oversampling_factor', 25, ...      % Oversampling factor
            'tx_gain', 0, ...                   % TX Gain [dB] (Range: -89.75 to 0 dB)
            'rx_gain', 50, ...                  % RX Gain [dB] (Range: -4 to 71 dB)
            'cable_length', 12, ...             % Cable length [m]
            'cable_loss_per_meter', 0.1, ...    % Cable loss [dB/m]
            'tx_power', 7, ...                  % Output power [dBm] (Range: -10 to +7 dBm)
            'antenna_diameter', 1.2, ...        % Parabolic antenna diameter [m]
            'antenna_efficiency', 0.6 ...       % Antenna efficiency
            );

        % Transponder and LO parameters (only used if useSameFrequency = false)
        transponderShift = 8089.5e6; % [Hz]
        Lo = 9750e6;                 % [Hz]

        % Flag: If true, RX center frequency is set equal to TX center frequency
        useSameFrequency = false;

        % Satellite parameters (assumptions)
        G_sat_rx = 30;    % Satellite receive antenna gain (Uplink) [dB]
        sat_rx_loss = 2;  % Satellite reception losses [dB]
        sat_amp_gain = 30;% Transponder amplification gain [dB]
        sat_tx_loss = 2;  % Satellite transmission losses [dB]
        G_sat_tx = 30;    % Satellite transmit antenna gain (Downlink) [dB]

        % Link parameters
        R = 35786e3;      % Distance to satellite [m]

        % Rain attenuation parameters (simplified ITU-R P.618 model)
        R_rain = 25;      % Rain rate in mm/h
        El_deg = 30;      % Elevation angle in degrees
        h_r = 5000;       % Effective rain height in m
        k_r = 0.0101;     % ITU-R P.838 coefficient
        alpha_r = 1.276;  % ITU-R P.838 coefficient
        reduction_factor = 0.6; % Reduction factor for effective rain path length

        % Receiver parameters
        T = 290;          % System temperature in Kelvin
        B = 1e6;          % Receiver bandwidth in Hz
    end

    properties (Dependent)
        lambda_u           % Wavelength for uplink [m]
        rx_center_frequency % Downlink center frequency [Hz]
        lambda_d           % Wavelength for downlink [m]
        cable_loss         % Total cable loss [dB]
        G_ground_u         % Base station antenna gain (Uplink) [dBi]
        G_ground_d         % Base station antenna gain (Downlink) [dBi]
        FSPL_u             % Free space path loss (Uplink) [dB]
        FSPL_d             % Free space path loss (Downlink) [dB]
        A_rain             % Rain attenuation [dB]
        EIRP_ground_u      % EIRP of base station (Uplink) [dBm]
        P_sat_rx           % Received power at satellite (Uplink) [dBm]
        P_sat_tx           % Output power of satellite (after transponder) [dBm]
        EIRP_sat           % EIRP of satellite (Downlink) [dBm]
        P_rx_ground        % Received power at base station (Downlink) [dBm]
        N_dBm              % Receiver noise power [dBm]
        SNR_dB             % Signal-to-noise ratio in dB
        SNR_linear         % Signal-to-noise ratio (linear)
    end

    methods
        %% Constructor: Overwrites parameters via name-value pairs
        function obj = SatelliteLink(varargin)
            if nargin > 0
                for k = 1:2:length(varargin)
                    if isfield(obj.baseStation, varargin{k})
                        obj.baseStation.(varargin{k}) = varargin{k+1};
                    elseif isprop(obj, varargin{k})
                        obj.(varargin{k}) = varargin{k+1};
                    end
                end
            end
        end

        %% Getters for dependent properties

        % Uplink wavelength (based on TX center frequency)
        function lambda = get.lambda_u(obj)
            lambda = 3e8 / obj.baseStation.tx_center_frequency;
        end

        % RX center frequency: either equal to TX center frequency or
        % calculated as tx_center_frequency + transponderShift - Lo
        function f_rx = get.rx_center_frequency(obj)
            if obj.useSameFrequency
                f_rx = obj.baseStation.tx_center_frequency;
            else
                f_rx = obj.baseStation.tx_center_frequency + obj.transponderShift - obj.Lo;
            end
        end

        % Downlink wavelength
        function lambda = get.lambda_d(obj)
            lambda = 3e8 / obj.rx_center_frequency;
        end

        % Total cable loss in dB
        function loss = get.cable_loss(obj)
            loss = obj.baseStation.cable_length * obj.baseStation.cable_loss_per_meter;
        end

        % Base station antenna gain for uplink
        function G = get.G_ground_u(obj)
            D = obj.baseStation.antenna_diameter;
            eta = obj.baseStation.antenna_efficiency;
            G = 10*log10(eta * (pi * D / obj.lambda_u)^2);
        end

        % Base station antenna gain for downlink
        function G = get.G_ground_d(obj)
            D = obj.baseStation.antenna_diameter;
            eta = obj.baseStation.antenna_efficiency;
            G = 10*log10(eta * (pi * D / obj.lambda_d)^2);
        end

        % Free space path loss (Uplink)
        function fspl = get.FSPL_u(obj)
            fspl = 20*log10(4*pi*obj.R/obj.lambda_u);
        end

        % Free space path loss (Downlink)
        function fspl = get.FSPL_d(obj)
            fspl = 20*log10(4*pi*obj.R/obj.lambda_d);
        end

        % Rain attenuation according to simplified ITU-R P.618 model
        function A = get.A_rain(obj)
            El = deg2rad(obj.El_deg);
            L_r = obj.h_r / sin(El);    % Geometric rain path length [m]
            L_r_km = L_r / 1000;        % in km
            gamma_R = obj.k_r * (obj.R_rain)^obj.alpha_r;  % dB/km
            L_eff = L_r_km * obj.reduction_factor;         % effective rain path length in km
            A = gamma_R * L_eff;
        end

        % EIRP of the base station for uplink
        function EIRP = get.EIRP_ground_u(obj)
            % EIRP = TX Power + TX Gain - Cable loss + Antenna gain
            EIRP = obj.baseStation.tx_power + obj.baseStation.tx_gain - obj.cable_loss + obj.G_ground_u;
        end

        % Received power at the satellite (Uplink)
        function P_rx = get.P_sat_rx(obj)
            P_rx = obj.EIRP_ground_u + obj.G_sat_rx - obj.FSPL_u - obj.A_rain - obj.sat_rx_loss;
        end

        % Output power of the satellite after transponder
        function P_tx = get.P_sat_tx(obj)
            P_tx = obj.P_sat_rx + obj.sat_amp_gain - obj.sat_tx_loss;
        end

        % EIRP of the satellite (Downlink)
        function EIRP_sat_val = get.EIRP_sat(obj)
            EIRP_sat_val = obj.P_sat_tx + obj.G_sat_tx;
        end

        % Received power at the base station (Downlink)
        function P_rx = get.P_rx_ground(obj)
            P_rx = obj.EIRP_sat + obj.G_ground_d - obj.cable_loss - obj.FSPL_d - obj.A_rain;
        end

        % Receiver noise power in dBm
        function N = get.N_dBm(obj)
            N_W = 1.38e-23 * obj.T * obj.B;
            N = 10*log10(N_W) + 30;
        end

        % Signal-to-noise ratio in dB
        function snr = get.SNR_dB(obj)
            snr = obj.P_rx_ground - obj.N_dBm;
        end

        % Signal-to-noise ratio (linear)
        function snr_lin = get.SNR_linear(obj)
            snr_lin = 10^(obj.SNR_dB/10);
        end

        %% Output of link budget results
        function displayLinkBudget(obj)
            fprintf('Link Budget Results\n');
            fprintf('  Base Station TX Center Frequency: %.2f Hz\n', obj.baseStation.tx_center_frequency);
            fprintf('  Calculated RX (Downlink) Center Frequency: %.2f Hz\n', obj.rx_center_frequency);
            fprintf('  Baseband Sample Rate: %.2f Hz\n', obj.baseStation.baseband_sample_rate);
            fprintf('  TX Gain: %.2f dB\n', obj.baseStation.tx_gain);
            fprintf('  RX Gain: %.2f dB\n', obj.baseStation.rx_gain);
            fprintf('  Antenna Gain (Uplink): %.2f dBi\n', obj.G_ground_u);
            fprintf('  Antenna Gain (Downlink): %.2f dBi\n', obj.G_ground_d);
            fprintf('  Cable Loss: %.2f dB\n', obj.cable_loss);
            fprintf('  FSPL Uplink: %.2f dB\n', obj.FSPL_u);
            fprintf('  FSPL Downlink: %.2f dB\n', obj.FSPL_d);
            fprintf('  Rain Attenuation: %.2f dB\n', obj.A_rain);
            fprintf('  EIRP Base Station (Uplink): %.2f dBm\n', obj.EIRP_ground_u);
            fprintf('  Power Received at Satellite (Uplink): %.2f dBm\n', obj.P_sat_rx);
            fprintf('  EIRP Satellite (Downlink): %.2f dBm\n', obj.EIRP_sat);
            fprintf('  Power Received at Base Station (Downlink): %.2f dBm\n', obj.P_rx_ground);
            fprintf('  Receiver Noise Power: %.2f dBm\n', obj.N_dBm);
            fprintf('  SNR: %.2f dB (%.2e linear)\n', obj.SNR_dB, obj.SNR_linear);
        end
    end
end
