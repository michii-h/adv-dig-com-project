function [vfcReceiveSignal] = simulate_qo100_channel(vfcTransmitSignal, stSat, stOFDM)
    % SIMULATE_QO100_CHANNEL Simulates the QO-100 satellite channel
    %   This function applies realistic satellite channel effects to a transmitted signal
    %   using p681LMSChannel model from the Satellite Communications Toolbox
    %
    %   Parameters:
    %   vfcTransmitSignal - Complex baseband signal to be transmitted
    %   stSat - Structure with satellite parameters from init_qo100_params
    %   stOFDM - Structure with OFDM parameters
    %
    %   Returns:
    %   vfcReceiveSignal - Complex baseband signal after passing through channel

    % Check if Satellite Communications Toolbox is available
    if ~license('test', 'Satellite_Comm_Toolbox')
        error('Satellite Communications Toolbox is required for this simulation');
    end

    % Get sampling rate from OFDM parameters or use default
    if isfield(stOFDM, 'fs')
        fs = stOFDM.fs;
    else
        fs = 8000; % Default sampling rate if not provided
        warning('Sampling rate not provided in stOFDM structure. Using default: 8000 Hz');
    end

    % Create a p681LMSChannel object for the uplink (2.4 GHz)
    uplink_channel = p681LMSChannel;
    uplink_channel.SampleRate = fs;
    uplink_channel.CarrierFrequency = stSat.uplink_freq;
    uplink_channel.ElevationAngle = stSat.el;            % Elevation angle from ground station
    uplink_channel.Environment = 'Suburban';         % Typical environment type

    % Create a p681LMSChannel object for the downlink (10.5 GHz)
    downlink_channel = p681LMSChannel;
    downlink_channel.SampleRate = fs;
    downlink_channel.CarrierFrequency = stSat.downlink_freq;
    downlink_channel.ElevationAngle = stSat.el;          % Same elevation angle
    downlink_channel.Environment = 'Suburban';       % Typical environment type

    % 1. Apply uplink effects

    % 1.1 Calculate and apply path loss for uplink
    uplink_attenuation = 10^(-stSat.path_loss_up/20);    % Convert dB to linear
    uplink_signal = vfcTransmitSignal * uplink_attenuation;

    % 1.2 Apply the p681LMSChannel model for uplink
    [uplink_channel_signal, uplink_path_gain] = uplink_channel(uplink_signal);

    % 2. Apply satellite transponder effects

    % 2.1 Apply nonlinear satellite transponder model
    transponder_signal = applyTransponderNonlinearity(uplink_channel_signal, stSat);

    % 2.2 Apply transponder gain
    transponder_gain = 10^(stSat.sat_transponder_gain/20);  % Convert dB to linear
    transponder_signal = transponder_signal * transponder_gain;

    % 3. Apply downlink effects

    % 3.1 Calculate and apply path loss for downlink
    downlink_attenuation = 10^(-stSat.path_loss_down/20);  % Convert dB to linear
    downlink_signal = transponder_signal * downlink_attenuation;

    % 3.2 Apply the p681LMSChannel model for downlink
    [downlink_channel_signal, downlink_path_gain] = downlink_channel(downlink_signal);

    % 4. Add thermal noise based on CNR
    signal_power = mean(abs(downlink_channel_signal).^2);
    CNR_linear = 10^(stSat.CNR/10);
    noise_power = signal_power / CNR_linear;

    % Generate complex Gaussian noise
    signal_length = length(downlink_channel_signal);
    noise = sqrt(noise_power/2) * (randn(signal_length, 1) + 1j*randn(signal_length, 1));

    % 5. Apply frequency offset due to oscillator drift
    % (Doppler effects are already included in the p681LMSChannel model)
    freq_offset_Hz = calculateOscillatorOffset(stSat);
    t = (0:signal_length-1)' / fs;
    freq_offset_factor = exp(1j * 2 * pi * freq_offset_Hz * t);

    % Apply frequency offset and add noise
    vfcReceiveSignal = downlink_channel_signal .* freq_offset_factor + noise;

    % Print channel information for debugging
    fprintf('QO-100 Channel Simulation:\n');
    fprintf('  Uplink Path Gain: %.2f dB\n', 20*log10(mean(abs(uplink_path_gain))));
    fprintf('  Downlink Path Gain: %.2f dB\n', 20*log10(mean(abs(downlink_path_gain))));
    fprintf('  Added Frequency Offset: %.2f Hz\n', freq_offset_Hz);
    fprintf('  Final CNR: %.2f dB\n', 10*log10(signal_power/noise_power));
end

function output_signal = applyTransponderNonlinearity(input_signal, stSat)
    % Apply nonlinear transponder characteristics (AM/AM and AM/PM conversion)

    % Get input signal envelope
    input_envelope = abs(input_signal);
    input_phase = angle(input_signal);

    % Normalize input power to transponder specs
    norm_factor = sqrt(mean(input_envelope.^2));
    norm_envelope = input_envelope / norm_factor;

    % Apply Saleh model for AM/AM conversion
    % Parameters can be adjusted based on actual transponder characteristics
    alpha_a = 2.1587; % AM/AM parameter
    beta_a = 1.1517;  % AM/AM parameter

    % AM/AM conversion
    output_magnitude = (alpha_a * norm_envelope) ./ (1 + beta_a * norm_envelope.^2);

    % Apply Saleh model for AM/PM conversion
    % Parameters can be adjusted based on actual transponder characteristics
    alpha_p = 4.0033; % AM/PM parameter
    beta_p = 9.1040;  % AM/PM parameter

    % AM/PM conversion (phase distortion based on input amplitude)
    phase_distortion = (alpha_p * norm_envelope.^2) ./ (1 + beta_p * norm_envelope.^2);
    output_phase = input_phase + phase_distortion;

    % Recombine magnitude and phase
    output_signal = output_magnitude .* exp(1j * output_phase);

    % Scale back to original power level
    output_signal = output_signal * norm_factor;
end

function freq_offset = calculateOscillatorOffset(stSat)
    % Calculate frequency offset due to oscillator instability
    % (Doppler is handled by the p681LMSChannel model)

    % Typical frequency stability values for oscillators
    uplink_oscillator_stability = 1e-6;  % 1 ppm
    downlink_oscillator_stability = 1e-6; % 1 ppm

    % Oscillator frequency errors
    uplink_osc_error = uplink_oscillator_stability * stSat.uplink_freq;
    downlink_osc_error = downlink_oscillator_stability * stSat.downlink_freq;

    % Total frequency offset from oscillator instabilities
    freq_offset = uplink_osc_error + downlink_osc_error;

    % Add some randomness to model variations
    freq_offset = freq_offset * (1 + 0.1 * (rand-0.5));
end
