function [vfcCaptureBuffer] = LoopbackAdalmPluto(vfcTransmitSignal, link)
    % Ensure Column Vector
    x = vfcTransmitSignal(:);

    % Use a conservative scale factor to avoid saturation
    powerScaleFactor = 0.9;
    txWaveform = x.*(1/max(abs(x))*powerScaleFactor);

    % Get parameters from link object
    bandwidth = link.baseStation.baseband_sample_rate;
    iOsf = link.baseStation.oversampling_factor;

    % Calculate required sample rate based on bandwidth
    % Using a factor of 2.5x bandwidth to ensure Nyquist and some margin
    requiredFs = bandwidth * 2.5;

    % Check if the sample rate is compatible with the bandwidth requirement
    if link.baseStation.baseband_sample_rate > requiredFs
        warning('Sample rate (%.1f kHz) may be too high for %.1f kHz bandwidth. Consider reducing to ~%.1f kHz.',...
                link.baseStation.baseband_sample_rate/1e3, bandwidth/1e3, requiredFs/1e3);
    end

    txWaveform = resample(txWaveform, iOsf, 1);

    % Init Transmission
    fs = link.baseStation.baseband_sample_rate*iOsf;     % Symbol rate with oversampling

    % Set TX and RX frequencies from the link object
    fcTx = link.baseStation.tx_center_frequency;
    fcRx = link.rx_center_frequency;

    % Get TX and RX gain from link object
    TxGain = link.baseStation.tx_gain;
    RxGain = link.baseStation.rx_gain;

    fprintf('\nUsing Adalm Pluto ...\n');

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
    fprintf('\nStarting transmission.\n')
    transmitRepeat(tx,txWaveform);

    captureLength = 2*length(x)*iOsf;

    fprintf('\nStarting capturing.\n')
    vfcCaptureBuffer = capture(rx, captureLength, 'Samples');

    vfcCaptureBuffer = resample(vfcCaptureBuffer,1,iOsf);

    release(tx);
    release(rx);
end
