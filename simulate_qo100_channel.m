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

    % Normalize input signal power to 1
    inputPower = sum(abs(vfcTransmitSignal).^2)/length(vfcTransmitSignal);
    vfcTransmitSignal = vfcTransmitSignal / sqrt(inputPower);

    % Apply the actual output power of the SDR based on our gain setting
    scaleFactor = sqrt(stSat.expected_output_power_w);
    vfcTransmitSignal = vfcTransmitSignal * scaleFactor;
    fprintf(' - Applied TX power scaling: %.2f dBm (%.6f W)\n', ...
            10*log10(stSat.expected_output_power_w)+30, stSat.expected_output_power_w);

    % Calculate free space path loss and apply it to the signal
    lambda = stSat.c / stSat.downlinkFreq;
    pathLoss = (4 * pi * stSat.slantRange / lambda)^2;
    pathLossDB = 10*log10(pathLoss);
    vfcSignal = vfcTransmitSignal / sqrt(pathLoss);
    fprintf(' - Applied free space path loss: %.2f dB\n', pathLossDB);

    % Calculate and apply the combined effects of antenna gain and other system gains/losses
    systemGainFactor = 10^(stSat.linkMargin/20); % Convert link margin to amplitude factor
    vfcSignal = vfcSignal * systemGainFactor;
    fprintf(' - Applied system gains/losses: %.2f dB\n', stSat.linkMargin);

    % Add thermal noise based on expected SNR
    signalPower = sum(abs(vfcSignal).^2)/length(vfcSignal);
    snr = 10^(stSat.expectedSNR/10); % Convert SNR from dB to linear
    noisePower = signalPower / snr;
    noise = sqrt(noisePower/2) * (randn(length(vfcSignal),1) + 1j*randn(length(vfcSignal),1));
    vfcSignal = vfcSignal + noise;
    fprintf(' - Added thermal noise for SNR: %.2f dB\n', stSat.expectedSNR);

    % 3. Apply phase noise
    % Use configurable phase noise or default to 0.01
    phaseNoiseVariance = 0.01; % Default value
    if isfield(stSat, 'phaseNoiseVariance')
        phaseNoiseVariance = stSat.phaseNoiseVariance;
    end
    phaseNoise = sqrt(phaseNoiseVariance) * randn(length(vfcSignal), 1);
    vfcSignal = vfcSignal .* exp(1j * phaseNoise);
    fprintf(' - Applied phase noise with variance: %.4f\n', phaseNoiseVariance);

    % 4. Apply frequency offset using proper time vector based on sample rate
    % QO-100 is geostationary so frequency offset is minimal but still present
    freqOffsetHz = 50; % Default value
    if isfield(stSat, 'freqOffset')
        freqOffsetHz = stSat.freqOffset;
    end
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
    fprintf(' - Applying nonlinear transponder effects (Saleh model)\n');

    % Default Saleh model parameters (can be overridden in stSat)
    alpha_a = 2.1587;
    beta_a = 1.1517;
    alpha_p = 4.0033;
    beta_p = 9.1040;

    % Allow parameter customization via stSat
    if isfield(stSat, 'saleh')
        if isfield(stSat.saleh, 'alpha_a'), alpha_a = stSat.saleh.alpha_a; end
        if isfield(stSat.saleh, 'beta_a'), beta_a = stSat.saleh.beta_a; end
        if isfield(stSat.saleh, 'alpha_p'), alpha_p = stSat.saleh.alpha_p; end
        if isfield(stSat.saleh, 'beta_p'), beta_p = stSat.saleh.beta_p; end
    end

    % Skip nonlinear modeling if specified in stSat
    if ~isfield(stSat, 'skipNonlinear') || ~stSat.skipNonlinear
        amplitude = abs(vfcSignal);
        phase = angle(vfcSignal);

        % Apply AM/AM conversion
        gainAM = alpha_a * amplitude ./ (1 + beta_a * amplitude.^2);
        % Apply AM/PM conversion (in radians)
        phasePM = alpha_p * amplitude.^2 ./ (1 + beta_p * amplitude.^2);

        % Apply the nonlinear distortion
        vfcSignal = gainAM .* exp(1j * (phase + phasePM));
        fprintf('   - Applied Saleh model with params: alpha_a=%.4f, beta_a=%.4f, alpha_p=%.4f, beta_p=%.4f\n', ...
            alpha_a, beta_a, alpha_p, beta_p);
    else
        fprintf('   - Nonlinear modeling skipped as requested\n');
    end

    % 6. Add transponder filtering
    if isfield(stSat, 'sampleRate') && length(vfcSignal) > 64
        if ~isfield(stSat, 'skipFiltering') || ~stSat.skipFiltering
            fprintf(' - Applying transponder filtering\n');
            actualSampleRate = stSat.sampleRate;
            if isfield(stSat, 'oversampling_factor')
                actualSampleRate = actualSampleRate * stSat.oversampling_factor;
            end

            % Get filter order, defaulting to 64 if not specified
            filterOrder = 64;
            if isfield(stSat, 'filterOrder')
                filterOrder = stSat.filterOrder;
            end

            % Default to transponderBW or 2.7 kHz if not specified
            filterBW = 2.7e3;
            if isfield(stSat, 'transponderBW')
                filterBW = stSat.transponderBW;
            end

            normalizedCutoff = min(0.95, filterBW / (actualSampleRate/2));
            h = fir1(filterOrder, normalizedCutoff);
            vfcSignal = filter(h, 1, vfcSignal);
            fprintf('   - Applied FIR filter with order %d and cutoff %.4f\n', filterOrder, normalizedCutoff);
        else
            fprintf('   - Transponder filtering skipped as requested\n');
        end
    end

    % Make sure it's a column vector
    vfcReceiveSignal = vfcSignal(:);

    fprintf('QO-100 channel simulation completed.\n');
end
