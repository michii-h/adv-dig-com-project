function vfcReceiveSignal = channel_sim(vfcTransmitSignal,stChannel)
    % Apply sample shift if needed
    iSampleShift = 0;
    vfcTransmitSignal = circshift(vfcTransmitSignal,[iSampleShift 0]);

    % Frequency offset (apply continuous rotation over time)
    vfcPhaser = exp(1j * 2 * pi * stChannel.fFreqOffset * [0:length(vfcTransmitSignal)-1]' / length(vfcTransmitSignal));
    vfcReceiveSignal = vfcTransmitSignal .* vfcPhaser;

    % Phase offset (constant phase shift)
    vfcReceiveSignal = vfcReceiveSignal * exp(1j * stChannel.fPhaseOffset);

    % Sample offset (fractional timing error)
    vfcFFT = fft(vfcReceiveSignal);
    vOmega = 2*pi*[0:length(vfcFFT)-1].'/length(vfcFFT);
    vfcReceiveSignal = ifft(vfcFFT .* fftshift(exp(1j * vOmega * stChannel.fSampleOffset)));

    % AWGN Channel
    vfcReceiveSignal = awgn(vfcReceiveSignal, stChannel.fSNRdB);

    % Multipath channel
    if length(stChannel.vfImpulseResponse) > 1
        vfcReceiveSignal = conv(vfcReceiveSignal, stChannel.vfImpulseResponse, 'same');
    end
end
