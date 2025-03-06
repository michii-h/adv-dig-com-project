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
    % For Adalm Pluto: 0 dB is maximum power, -89.75 dB is minimum
    if stSat.adalm_txGain > 0
        warning('Adalm Pluto TX gain should be <= 0 dB. Using 0 dB (maximum power).');
        scaleFactor = 1.0;  % No attenuation at maximum power
    else
        % Convert negative dB to amplitude scaling (negative dB means attenuation)
        scaleFactor = 10^(stSat.adalm_txGain/20);
    end
    vfcTransmitSignal = vfcTransmitSignal * scaleFactor;
    fprintf(' - Applied TX gain scaling: %.2f dB (factor: %.4f)\n', stSat.adalm_txGain, scaleFactor);

    % 1. Apply path loss and power scaling based on link budget
    linkScaleFactor = 10^(-stSat.linkMargin/20); % Convert dB to amplitude scaling
    vfcSignal = vfcTransmitSignal * linkScaleFactor;
    fprintf(' - Applied link margin scaling: %.2f dB (factor: %.4f)\n', -stSat.linkMargin, linkScaleFactor);

    % 2. Add thermal noise based on SNR derived from link margin
    signalPower = sum(abs(vfcSignal).^2)/length(vfcSignal);
    snr = 10^(stSat.expectedSNR/10); % Convert SNR from dB to linear
    noisePower = signalPower / snr;
    noise = sqrt(noisePower/2) * (randn(length(vfcSignal),1) + 1j*randn(length(vfcSignal),1));
    vfcSignal = vfcSignal + noise;
    fprintf(' - Added thermal noise for SNR: %.2f dB\n', stSat.expectedSNR);

    % 3. Apply phase noise
    phaseNoiseVariance = 0.01; % Adjust based on satellite transponder characteristics
    phaseNoise = sqrt(phaseNoiseVariance) * randn(length(vfcSignal), 1);
    vfcSignal = vfcSignal .* exp(1j * phaseNoise);
    fprintf(' - Applied phase noise with variance: %.4f\n', phaseNoiseVariance);

    % 4. Apply frequency offset using proper time vector based on sample rate
    % QO-100 is geostationary so frequency offset is minimal but still present
    freqOffsetHz = 50; % Hz
    if isfield(stSat, 'sampleRate')
        actualSampleRate = stSat.sampleRate;
        if isfield(stSat, 'oversampling_factor')
            actualSampleRate = actualSampleRate * stSat.oversampling_factor;
        end
        t = (0:length(vfcSignal)-1)' / actualSampleRate;
        vfcSignal = vfcSignal .* exp(1j * 2 * pi * freqOffsetHz * t);
        fprintf(' - Applied frequency offset: %.2f Hz\n', freqOffsetHz);
    else
        warning('Sample rate not specified; using normalized frequency offset');
        freqOffset = 0.001; % Normalized frequency offset [0-1]
        t = (0:length(vfcSignal)-1)' / length(vfcSignal);
        vfcSignal = vfcSignal .* exp(1j * 2 * pi * freqOffset * t);
        fprintf(' - Applied normalized frequency offset: %.4f\n', freqOffset);
    end

    % 5. Model nonlinear transponder effects (Saleh model)
    % fprintf(' - Applying nonlinear transponder effects (Saleh model)\n');
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

    % 6. Add transponder filtering
    % if isfield(stSat, 'transponderBW') && length(vfcSignal) > 64
    %     fprintf(' - Applying transponder filtering\n');
    %     if isfield(stSat, 'sampleRate')
    %         actualSampleRate = stSat.sampleRate;
    %         if isfield(stSat, 'oversampling_factor')
    %             actualSampleRate = actualSampleRate * stSat.oversampling_factor;
    %         end
    %         normalizedCutoff = min(0.95, stSat.transponderBW / (actualSampleRate/2));
    %         h = fir1(64, normalizedCutoff);
    %         vfcSignal = filter(h, 1, vfcSignal);
    %     end
    % end

    % Make sure it's a column vector
    vfcReceiveSignal = vfcSignal(:);

    fprintf('QO-100 channel simulation completed.\n');
end
