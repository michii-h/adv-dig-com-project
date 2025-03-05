function vfcReceiveSignal = simulate_qo100_channel(vfcTransmitSignal, stSat)
    % SIMULATE_QO100_CHANNEL Simulates the QO-100 satellite communication channel
    % This function uses the MATLAB Satellite Communications Toolbox to model
    % the channel effects for the Qatar-OSCAR 100 satellite.
    %
    % Inputs:
    %   vfcTransmitSignal - Complex baseband transmit signal
    %   stSat - Structure containing QO-100 parameters

    % Ensure column vector
    vfcTransmitSignal = vfcTransmitSignal(:);

    fprintf('Simulating QO-100 satellite channel with Link Margin: %.2f dB\n', stSat.linkMargin);

    % 0. Apply tx power scaling based on Adalm Pluto gain
    scaleFactor = 10^(-stSat.adalm_txGain/20); % Convert dB to amplitude scaling; Adalm Pluto gain is negative
    vfcTransmitSignal = vfcTransmitSignal * scaleFactor;

    % 1. Apply path loss and power scaling based on link budget
    scaleFactor = 10^(-stSat.linkMargin/20); % Convert dB to amplitude scaling
    vfcSignal = vfcTransmitSignal * scaleFactor;

    % 2. Add thermal noise based on SNR derived from link margin
    signalPower = sum(abs(vfcSignal).^2)/length(vfcSignal);
    snr = 10^(stSat.expectedSNR/10); % Convert SNR from dB to linear
    noisePower = signalPower / snr;
    noise = sqrt(noisePower/2) * (randn(length(vfcSignal),1) + 1j*randn(length(vfcSignal),1));
    vfcSignal = vfcSignal + noise;

    % 3. Apply phase noise
    phaseNoiseVariance = 0.01; % Adjust based on satellite transponder characteristics
    phaseNoise = sqrt(phaseNoiseVariance) * randn(length(vfcSignal), 1);
    vfcSignal = vfcSignal .* exp(1j * phaseNoise);

    % 4. Apply frequency offset (minimal for geostationary satellites)
    % freqOffset = 50; % Hz, typical frequency offset
    % t = (0:length(vfcSignal)-1)' / stSat.sampleRate;
    freqOffset = 0.001; % Frequency offset [fs/iNfft]
    t = ( 0:length(vfcTransmitSignal)-1 )' / length(vfcTransmitSignal);
    vfcSignal = vfcSignal .* exp(1j * 2 * pi * freqOffset * t);

    % % 5. Model nonlinear transponder effects (simplified Saleh model)
    % % AM/AM distortion
    % alpha_a = 2.1587;
    % beta_a = 1.1517;
    % % AM/PM distortion
    % alpha_p = 4.0033;
    % beta_p = 9.1040;

    % amplitude = abs(vfcSignal);
    % phase = angle(vfcSignal);

    % % Apply AM/AM conversion
    % gainAM = alpha_a * amplitude ./ (1 + beta_a * amplitude.^2);
    % % Apply AM/PM conversion (in radians)
    % phasePM = alpha_p * amplitude.^2 ./ (1 + beta_p * amplitude.^2);

    % % Apply the nonlinear distortion
    % vfcSignal = gainAM .* exp(1j * (phase + phasePM));

    % % 6. Add transponder filtering (simplified)
    % if length(vfcSignal) > 64
    %     h = fir1(64, stSat.transponderBW/stSat.sampleRate);
    %     vfcSignal = filter(h, 1, vfcSignal);
    % end

    % Make sure it's a column vector
    vfcReceiveSignal = vfcSignal(:);
end
