function [vfcCaptureBuffer] = LoopbackAdalmPluto(vfcTransmitSignal,stZedboard)
  % Scale the normalized signal to avoid saturation of RF stages
        x = vfcTransmitSignal(:);        
        powerScaleFactor = 0.8;
        txWaveform = x.*(1/max(abs(x))*powerScaleFactor);

        iOsf = 50;

        txWaveform = resample(txWaveform,iOsf,1);

        
        % Init Transmission
        
        fs = stZedboard.fs*iOsf;     % Symbol rate
        fc = stZedboard.fc; % Carrier center Frequency
        TxGain = stZedboard.TxGain; % Gain, specified as a scalar from -89.75 to 0 dB with a resolution of 0.25 dB
        RxGain = stZedboard.RxGain;  % Radio receiver gain in dB, specified as a scalar from -4 to 71
        chIPAddress = stZedboard.IPAddress;
        vTxChannelMapping = stZedboard.TxChannelMapping;
        vRxChannelMapping = stZedboard.RxChannelMapping;

        % Set up TX Radio    
        %tx = sdrtx('FMCOMMS5', 'IPAddress',chIPAddress);
        tx = sdrtx('AD936x', 'IPAddress',chIPAddress);
    
    
        
        tx.CenterFrequency      = fc;
        tx.BasebandSampleRate   = fs;
        tx.Gain                 = TxGain;
        tx.ChannelMapping  = vTxChannelMapping;
        
        % Set up RX Radio
        %rx = sdrrx('FMCOMMS5', 'IPAddress',chIPAddress);
        rx = sdrrx('AD936x', 'IPAddress',chIPAddress);
        rx.BasebandSampleRate   = tx.BasebandSampleRate;
        rx.CenterFrequency      = tx.CenterFrequency;
        rx.GainSource           = 'Manual';
        rx.Gain                 = RxGain;
        rx.OutputDataType       = 'double';
        rx.ChannelMapping  = vRxChannelMapping;

        % Pass data through radio
        fprintf('\nStarting transmission.\n')
        vTxWaveform=repmat(txWaveform,1,size(vTxChannelMapping,2));
        


        transmitRepeat(tx,vTxWaveform);
        
        captureLength = 2*length(x)*iOsf;
        fprintf('\nStarting capturing.\n')

        vfcCaptureBuffer = capture(rx, captureLength, 'Samples');

        
        vfcCaptureBuffer = resample(vfcCaptureBuffer,1,iOsf);
        
        release(tx);
        release(rx);