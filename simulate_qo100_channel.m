function vfcReceiveSignal = simulate_qo100_channel(vfcTransmitSignal, link)
    % Simulates the QO-100 satellite communication channel
    %
    % Inputs:
    %   vfcTransmitSignal - Complex baseband transmit signal
    %   link - SatelliteLink object containing parameters

    % Ensure column vector
    vfcTransmitSignal = vfcTransmitSignal(:);

    fprintf('Simulating QO-100 satellite channel with SNR: %.2f dB\n', link.SNR_dB);

    % Normalize input signal power to 1
    % inputPower = sum(abs(vfcTransmitSignal).^2)/length(vfcTransmitSignal);
    % vfcTransmitSignal = vfcTransmitSignal / sqrt(inputPower);

    % Instead of applying all individual gains and losses, we can
    % directly use the calculated SNR from the link budget
    % The normalized signal already has a power of 1
    % So we can directly add noise based on the target SNR
    snr = 10^(link.SNR_dB/10); % Convert SNR from dB to linear
    noisePower = 1 / snr;
    noise = sqrt(noisePower/2) * (randn(length(vfcTransmitSignal),1) + 1j*randn(length(vfcTransmitSignal),1));
    vfcSignal = vfcTransmitSignal + noise;
    fprintf('  Added thermal noise for SNR: %.2f dB\n', link.SNR_dB);

    % 3. Apply phase noise
    % Use a default phase noise variance
    phaseNoiseVariance = 0.01;
    phaseNoise = sqrt(phaseNoiseVariance) * randn(length(vfcSignal), 1);
    vfcSignal = vfcSignal .* exp(1j * phaseNoise);
    fprintf('  Applied phase noise with variance: %.4f\n', phaseNoiseVariance);

    % 4. Apply frequency offset using proper time vector based on sample rate
    freqOffsetHz = 50; % Default frequency offset
    actualSampleRate = link.baseStation.baseband_sample_rate * link.baseStation.oversampling_factor;
    t = (0:length(vfcSignal)-1)' / actualSampleRate;
    vfcSignal = vfcSignal .* exp(1j * 2 * pi * freqOffsetHz * t);
    fprintf('  Applied frequency offset: %.2f Hz\n', freqOffsetHz);

    % Make sure it's a column vector
    vfcReceiveSignal = vfcSignal(:);

    fprintf('QO-100 channel simulation completed.\n');
end
