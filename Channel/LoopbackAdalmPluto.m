function [vfcCaptureBuffer] = LoopbackAdalmPluto(vfcTransmitSignal,stAdalmPluto)
  % Scale the normalized signal to avoid saturation of RF stages
        x = vfcTransmitSignal(:);        
        powerScaleFactor = 0.8;
        txWaveform = x.*(1/max(abs(x))*powerScaleFactor);

        iOsf = 10;

        txWaveform = resample(txWaveform,iOsf,1);

        
        % Init Transmission
        
        fs = stAdalmPluto.fs*iOsf;     % Symbol rate
        fc = stAdalmPluto.fc; % Carrier center Frequency
        TxGain = stAdalmPluto.TxGain; % Gain, specified as a scalar from -89.75 to 0 dB with a resolution of 0.25 dB
        RxGain = stAdalmPluto.RxGain;  % Radio receiver gain in dB, specified as a scalar from -4 to 71
        
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
        fprintf('\nStarting transmission.\n')
        transmitRepeat(tx,txWaveform);
        
        captureLength = 2*length(x)*iOsf;
        fprintf('\nStarting capturing.\n')

        vfcCaptureBuffer = capture(rx, captureLength, 'Samples');

        
        vfcCaptureBuffer = resample(vfcCaptureBuffer,1,iOsf);
        
        release(tx);
        release(rx);