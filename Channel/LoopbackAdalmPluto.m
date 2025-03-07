function [vfcCaptureBuffer] = LoopbackAdalmPluto(vfcTransmitSignal,stAdalmPluto, i)
    % Normalize input signal to avoid clipping
    x = vfcTransmitSignal(:);

    % Normalize input signal power to 1
    % inputPower = sum(abs(x).^2)/length(x);
    % x = x / sqrt(inputPower);

    % Use a more conservative scale factor to avoid saturation
    % The actual power will be determined by the TX gain setting
    powerScaleFactor = 0.8;
    txWaveform = x.*(1/max(abs(x))*powerScaleFactor);

    % Use oversampling factor from stAdalmPluto if available, otherwise default to 10
    if isfield(stAdalmPluto, 'oversampling_factor')
        iOsf = stAdalmPluto.oversampling_factor;
    else
        iOsf = 10;
    end

    txWaveform = resample(txWaveform,iOsf,1);

    % Init Transmission
    fs = stAdalmPluto.fs*iOsf;     % Symbol rate
    fc = stAdalmPluto.fc;          % Carrier center Frequency

    % Ensure TX gain is within valid range for Adalm Pluto (-89.75 to 0 dB)
    TxGain = min(0, max(-89.75, stAdalmPluto.TxGain));
    RxGain = min(71, max(-4, stAdalmPluto.RxGain));

    % Calculate approximate output power based on gain setting
    approxPowerDbm = stAdalmPluto.max_power_dbm + TxGain;  % Reduce by gain amount (negative dB)

    if i==1
        fprintf('Using Adalm Pluto with:\n');
        fprintf(' - Carrier: %.3f Hz\n', fc);
        fprintf(' - Sample Rate: %.3f MHz\n', fs/1e6);
        fprintf(' - TX Gain: %.2f dB (range: -89.75 to 0 dB, where 0 dB is max power)\n', TxGain);
        fprintf(' - Approximate TX Power: %.2f dBm\n', approxPowerDbm);
        fprintf(' - RX Gain: %.2f dB\n', RxGain);
    end

    % Set up TX Radio
    tx = sdrtx('Pluto');

    tx.CenterFrequency      = fc;
    tx.BasebandSampleRate   = fs;
    tx.Gain                 = TxGain;

    % Set up RX Radio
    rx = sdrrx('Pluto');
    rx.BasebandSampleRate   = tx.BasebandSampleRate;
    rx.CenterFrequency      = tx.CenterFrequency;
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
