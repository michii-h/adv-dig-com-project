function vfcReceiveSignal = simulate_qo100_channel(vfcTransmitSignal, link)
    % SIMULATE_QO100_CHANNEL Simulates the QO-100 satellite communication channel
    % This function uses the MATLAB Satellite Communications Toolbox to model
    % the channel effects for the Qatar-OSCAR 100 satellite.
    %
    % Inputs:
    %   vfcTransmitSignal - Complex baseband transmit signal
    %   link - SatelliteLink object containing parameters

    % Ensure column vector
    vfcTransmitSignal = vfcTransmitSignal(:);

    fprintf('Simulating QO-100 satellite channel with SNR: %.2f dB\n', link.SNR_dB);

    % Normalize input signal power to 1
    inputPower = sum(abs(vfcTransmitSignal).^2)/length(vfcTransmitSignal);
    vfcTransmitSignal = vfcTransmitSignal / sqrt(inputPower);

    % Apply the actual output power based on link settings
    % Convert from dBm to Watts
    outputPowerW = 10^((link.baseStation.tx_power - 30)/10);
    scaleFactor = sqrt(outputPowerW);
    vfcTransmitSignal = vfcTransmitSignal * scaleFactor;
    fprintf(' - Applied TX power scaling: %.2f dBm (%.6f W)\n', link.baseStation.tx_power, outputPowerW);

    % Calculate free space path loss and apply it to the signal
    lambda = 3e8 / link.baseStation.tx_center_frequency;
    pathLoss = (4 * pi * link.R / lambda)^2;
    pathLossDB = 10*log10(pathLoss);
    vfcSignal = vfcTransmitSignal / sqrt(pathLoss);
    fprintf(' - Applied free space path loss: %.2f dB\n', pathLossDB);

    % Calculate and apply the combined effects of antenna gain and other system gains/losses
    % Use SNR margin as a proxy for link margin
    systemGainFactor = 10^((link.SNR_dB - 20)/20); % Converting to amplitude factor
    vfcSignal = vfcSignal * systemGainFactor;
    fprintf(' - Applied system gains/losses based on SNR: %.2f dB\n', link.SNR_dB);

    % Add thermal noise based on expected SNR
    signalPower = sum(abs(vfcSignal).^2)/length(vfcSignal);
    snr = 10^(link.SNR_dB/10); % Convert SNR from dB to linear
    noisePower = signalPower / snr;
    noise = sqrt(noisePower/2) * (randn(length(vfcSignal),1) + 1j*randn(length(vfcSignal),1));
    vfcSignal = vfcSignal + noise;
    fprintf(' - Added thermal noise for SNR: %.2f dB\n', link.SNR_dB);

    % 3. Apply phase noise
    % Use a default phase noise variance
    phaseNoiseVariance = 0.01;
    phaseNoise = sqrt(phaseNoiseVariance) * randn(length(vfcSignal), 1);
    vfcSignal = vfcSignal .* exp(1j * phaseNoise);
    fprintf(' - Applied phase noise with variance: %.4f\n', phaseNoiseVariance);

    % 4. Apply frequency offset using proper time vector based on sample rate
    freqOffsetHz = 50; % Default frequency offset
    actualSampleRate = link.baseStation.baseband_sample_rate * link.baseStation.oversampling_factor;
    t = (0:length(vfcSignal)-1)' / actualSampleRate;
    vfcSignal = vfcSignal .* exp(1j * 2 * pi * freqOffsetHz * t);
    fprintf(' - Applied frequency offset: %.2f Hz\n', freqOffsetHz);

    % 5. Model nonlinear transponder effects (Saleh model)
    fprintf(' - Applying nonlinear transponder effects (Saleh model)\n');

    % Default Saleh model parameters
    alpha_a = 2.1587;
    beta_a = 1.1517;
    alpha_p = 4.0033;
    beta_p = 9.1040;

    % Apply the Saleh model
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

    % 6. Add transponder filtering
    fprintf(' - Applying transponder filtering\n');
    actualSampleRate = link.baseStation.baseband_sample_rate * link.baseStation.oversampling_factor;
    filterOrder = 64; % Default filter order
    filterBW = link.baseStation.baseband_sample_rate; % Use baseband sample rate as bandwidth
    normalizedCutoff = min(0.95, filterBW / (actualSampleRate/2));
    h = fir1(filterOrder, normalizedCutoff);
    vfcSignal = filter(h, 1, vfcSignal);
    fprintf('   - Applied FIR filter with order %d and cutoff %.4f\n', filterOrder, normalizedCutoff);

    % Make sure it's a column vector
    vfcReceiveSignal = vfcSignal(:);

    fprintf('QO-100 channel simulation completed.\n');
end
