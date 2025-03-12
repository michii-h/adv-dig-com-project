% 'tx_center_frequency', 2400.172 * 1e6, ...   % Uplink-Centerfrequenz [Hz]
% 'baseband_sample_rate', 2.7e3, ...   % Baseband-Samplerate [Hz]
% 'tx_gain', 0, ...                    % TX Gain [dB]
% 'rx_gain', 50, ...                   % RX Gain [dB]

classdef SatelliteLink
    % Diese Klasse berechnet das Link Budget für einen Satellitenlink unter
    % Berücksichtigung von Basisstationsparametern und Satellitenparametern.
    % Die Basisstationsparameter (TX-/RX‑Centerfrequenz, Baseband-Samplerate,
    % TX-/RX‑Gain, Ausgangsleistung, Kabelverluste und Antennendaten) sind in
    % einem eigenen Strukturfeld 'baseStation' zusammengefasst.
    %
    % Standardmäßig wird die RX‑Centerfrequenz berechnet als:
    %    rx_center_frequency = tx_center_frequency + transponderShift - Lo
    % Falls useSameFrequency = true gesetzt ist, wird stattdessen
    %    rx_center_frequency = tx_center_frequency
    % verwendet.

    properties
        % Basisstationsparameter als Struktur
        baseStation = struct( ...
            'tx_center_frequency', 2400.172 * 1e6, ...   % Uplink-Centerfrequenz [Hz]
            'baseband_sample_rate', 2.7e3, ...  % Baseband-Samplerate [Hz]
            'oversampling_factor', 25, ...      % Oversampling-Faktor
            'tx_gain', 0, ...                   % TX Gain [dB] (Bereich: -89.75 bis 0 dB)
            'rx_gain', 50, ...                  % RX Gain [dB] (Bereich: -4 bis 71 dB)
            'cable_length', 12, ...             % Kabellänge [m]
            'cable_loss_per_meter', 0.1, ...    % Kabelverlust [dB/m]
            'tx_power', 7, ...                  % Ausgangsleistung [dBm] (Bereich: -10 bis +7 dBm)
            'antenna_diameter', 1.2, ...        % Durchmesser der Parabolantenne [m]
            'antenna_efficiency', 0.6 ...       % Antenneneffizienz
            );

        % Transponder- und LO-Parameter (nur genutzt, falls useSameFrequency = false)
        transponderShift = 8089.5e6; % [Hz]
        Lo = 9750e6;                 % [Hz]

        % Flag: Wenn true, wird RX-Centerfrequenz gleich TX-Centerfrequenz gesetzt
        useSameFrequency = false;

        % Satellitenparameter (Annahmen)
        G_sat_rx = 30;    % Empfangsantennengewinn des Satelliten (Uplink) [dB]
        sat_rx_loss = 2;  % Empfangsverluste im Satelliten [dB]
        sat_amp_gain = 30;% Transponderverstärkung im Satelliten [dB]
        sat_tx_loss = 2;  % Sendeverluste im Satelliten [dB]
        G_sat_tx = 30;    % Satellit Sendantennengewinn (Downlink) [dB]

        % Linkparameter
        R = 35786e3;      % Entfernung zum Satelliten [m]

        % Regenabschwächungsparameter (vereinfachtes ITU‑R P.618 Modell)
        R_rain = 25;      % Regenrate in mm/h
        El_deg = 30;      % Elevationswinkel in Grad
        h_r = 5000;       % Effektive Regenhöhe in m
        k_r = 0.0101;     % ITU‑R P.838 Koeffizient
        alpha_r = 1.276;  % ITU‑R P.838 Koeffizient
        reduction_factor = 0.6; % Reduktionsfaktor für effektive Regenpfadlänge

        % Empfangsparameter
        T = 290;          % Systemtemperatur in Kelvin
        B = 1e6;          % Empfängerbandbreite in Hz
    end

    properties (Dependent)
        lambda_u           % Wellenlänge für Uplink [m]
        rx_center_frequency % Downlink-Centerfrequenz [Hz]
        lambda_d           % Wellenlänge für Downlink [m]
        cable_loss         % Gesamter Kabelverlust [dB]
        G_ground_u         % Antennengewinn der Basisstation (Uplink) [dBi]
        G_ground_d         % Antennengewinn der Basisstation (Downlink) [dBi]
        FSPL_u             % Freiraumdämpfung (Uplink) [dB]
        FSPL_d             % Freiraumdämpfung (Downlink) [dB]
        A_rain             % Regenabschwächung [dB]
        EIRP_ground_u      % EIRP der Basisstation (Uplink) [dBm]
        P_sat_rx           % Empfangene Leistung am Satelliten (Uplink) [dBm]
        P_sat_tx           % Ausgangsleistung des Satelliten (nach Transponder) [dBm]
        EIRP_sat           % EIRP des Satelliten (Downlink) [dBm]
        P_rx_ground        % Empfangene Leistung an der Basisstation (Downlink) [dBm]
        N_dBm              % Rauschleistung des Empfängers [dBm]
        SNR_dB             % Signal-Rausch-Verhältnis in dB
        SNR_linear         % Signal-Rausch-Verhältnis (linear)
    end

    methods
        %% Konstruktor: Überschreibt Parameter per Name-Value-Paar
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

        %% Getter für abhängige Eigenschaften

        % Uplink-Wellenlänge (basierend auf TX-Centerfrequenz)
        function lambda = get.lambda_u(obj)
            lambda = 3e8 / obj.baseStation.tx_center_frequency;
        end

        % RX-Centerfrequenz: entweder gleich TX-Centerfrequenz oder
        % berechnet als tx_center_frequency + transponderShift - Lo
        function f_rx = get.rx_center_frequency(obj)
            if obj.useSameFrequency
                f_rx = obj.baseStation.tx_center_frequency;
            else
                f_rx = obj.baseStation.tx_center_frequency + obj.transponderShift - obj.Lo;
            end
        end

        % Downlink-Wellenlänge
        function lambda = get.lambda_d(obj)
            lambda = 3e8 / obj.rx_center_frequency;
        end

        % Gesamter Kabelverlust in dB
        function loss = get.cable_loss(obj)
            loss = obj.baseStation.cable_length * obj.baseStation.cable_loss_per_meter;
        end

        % Antennengewinn der Basisstation im Uplink
        function G = get.G_ground_u(obj)
            D = obj.baseStation.antenna_diameter;
            eta = obj.baseStation.antenna_efficiency;
            G = 10*log10(eta * (pi * D / obj.lambda_u)^2);
        end

        % Antennengewinn der Basisstation im Downlink
        function G = get.G_ground_d(obj)
            D = obj.baseStation.antenna_diameter;
            eta = obj.baseStation.antenna_efficiency;
            G = 10*log10(eta * (pi * D / obj.lambda_d)^2);
        end

        % Freiraumdämpfung (Uplink)
        function fspl = get.FSPL_u(obj)
            fspl = 20*log10(4*pi*obj.R/obj.lambda_u);
        end

        % Freiraumdämpfung (Downlink)
        function fspl = get.FSPL_d(obj)
            fspl = 20*log10(4*pi*obj.R/obj.lambda_d);
        end

        % Regenabschwächung nach vereinfachtem ITU‑R P.618 Modell
        function A = get.A_rain(obj)
            El = deg2rad(obj.El_deg);
            L_r = obj.h_r / sin(El);    % Geometrische Regenpfadlänge [m]
            L_r_km = L_r / 1000;        % in km
            gamma_R = obj.k_r * (obj.R_rain)^obj.alpha_r;  % dB/km
            L_eff = L_r_km * obj.reduction_factor;         % effektive Regenpfadlänge in km
            A = gamma_R * L_eff;
        end

        % EIRP der Basisstation im Uplink
        function EIRP = get.EIRP_ground_u(obj)
            % EIRP = TX Power + TX Gain - Kabelverlust + Antennengewinn
            EIRP = obj.baseStation.tx_power + obj.baseStation.tx_gain - obj.cable_loss + obj.G_ground_u;
        end

        % Empfangene Leistung am Satelliten (Uplink)
        function P_rx = get.P_sat_rx(obj)
            P_rx = obj.EIRP_ground_u + obj.G_sat_rx - obj.FSPL_u - obj.A_rain - obj.sat_rx_loss;
        end

        % Ausgangsleistung des Satelliten nach Transponder
        function P_tx = get.P_sat_tx(obj)
            P_tx = obj.P_sat_rx + obj.sat_amp_gain - obj.sat_tx_loss;
        end

        % EIRP des Satelliten (Downlink)
        function EIRP_sat_val = get.EIRP_sat(obj)
            EIRP_sat_val = obj.P_sat_tx + obj.G_sat_tx;
        end

        % Empfangene Leistung an der Basisstation (Downlink)
        function P_rx = get.P_rx_ground(obj)
            P_rx = obj.EIRP_sat + obj.G_ground_d - obj.cable_loss - obj.FSPL_d - obj.A_rain;
        end

        % Rauschleistung des Empfängers in dBm
        function N = get.N_dBm(obj)
            N_W = 1.38e-23 * obj.T * obj.B;
            N = 10*log10(N_W) + 30;
        end

        % Signal-Rausch-Verhältnis in dB
        function snr = get.SNR_dB(obj)
            snr = obj.P_rx_ground - obj.N_dBm;
        end

        % Signal-Rausch-Verhältnis (linear)
        function snr_lin = get.SNR_linear(obj)
            snr_lin = 10^(obj.SNR_dB/10);
        end

        %% Ausgabe der Link Budget Ergebnisse
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
