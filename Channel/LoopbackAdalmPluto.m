function [vfcCaptureBuffer] = LoopbackAdalmPluto(vfcTransmitSignal, stSat, i)
    % Ensure Column Vector
    x = vfcTransmitSignal(:);

    % Use a conservative scale factor to avoid saturation
    powerScaleFactor = 0.9;
    txWaveform = x.*(1/max(abs(x))*powerScaleFactor);

    bandwidth = stSat.bandwidth;
    iOsf = stSat.oversampling_factor;

    % Calculate required sample rate based on bandwidth
    % Using a factor of 2.5x bandwidth to ensure Nyquist and some margin
    requiredFs = bandwidth * 2.5;

    % Check if the sample rate is compatible with the bandwidth requirement
    if stSat.fs > requiredFs
        warning('Sample rate (%.1f kHz) may be too high for %.1f kHz bandwidth. Consider reducing to ~%.1f kHz.',...
                stSat.fs/1e3, bandwidth/1e3, requiredFs/1e3);
    end

    txWaveform = resample(txWaveform, iOsf, 1);

    % Init Transmission
    fs = stSat.fs*iOsf;     % Symbol rate with oversampling

    % Set TX and RX frequencies
    if isfield(stSat, 'fcTx') && isfield(stSat, 'fcRx')
        fcTx = stSat.fcTx;
        fcRx = stSat.fcRx;
    else
        fcTx = stSat.fc;
        fcRx = stSat.fc;
    end

    % Ensure TX gain is within valid range for Adalm Pluto
    TxGain = stSat.txGain;
    RxGain = stSat.rxGain;

    % Calculate approximate output power based on gain setting
    approxPowerDbm = stSat.max_power_dbm + TxGain;  % Reduce by gain amount (negative dB)

    if i==1
        fprintf('Using Adalm Pluto for QO-100:\n');
        fprintf(' - TX Frequency: %.6f MHz\n', fcTx/1e6);
        fprintf(' - RX Frequency: %.6f MHz\n', fcRx/1e6);
        fprintf(' - Sample Rate: %.3f MHz\n', fs/1e6);
        fprintf(' - Bandwidth Limit: %.1f kHz\n', bandwidth/1e3);
        fprintf(' - TX Gain: %.2f dB (range: -89.75 to 0 dB, where 0 dB is max power)\n', TxGain);
        fprintf(' - Approximate TX Power: %.2f dBm\n', approxPowerDbm);
        fprintf(' - RX Gain: %.2f dB\n', RxGain);
    end

    % Set up TX Radio
    tx = sdrtx('Pluto');

    tx.CenterFrequency      = fcTx;
    tx.BasebandSampleRate   = fs;
    tx.Gain                 = TxGain;

    % Set up RX Radio
    rx = sdrrx('Pluto');
    rx.BasebandSampleRate   = tx.BasebandSampleRate;
    rx.CenterFrequency      = fcRx;
    rx.GainSource           = 'Manual';
    rx.Gain                 = RxGain;
    rx.OutputDataType       = 'double';

    % Pass data through radio
    if i==1
        fprintf('\nStarting transmission.\n')
    end
    transmitRepeat(tx,txWaveform);

    captureLength = 2*length(x)*iOsf;
    if i==1
        fprintf('\nStarting capturing.\n')
    end
    vfcCaptureBuffer = capture(rx, captureLength, 'Samples');

    vfcCaptureBuffer = resample(vfcCaptureBuffer,1,iOsf);

    release(tx);
    release(rx);
end
