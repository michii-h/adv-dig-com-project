function vfcReceiveSignal = channel_sim(vfcTransmitSignal,stChannel)
    % Simulate iSampleShift = 17 between Tx and Rx
    iSampleShift = 15;
    vfcTransmitSignal = circshift(vfcTransmitSignal,[iSampleShift 0]);


    % Frequency- and Phase-Offset
    vfcPhaser = exp(1j*stChannel.fFreqOffset*[0:length(vfcTransmitSignal)-1]+1j*stChannel.fPhaseOffset);
    vfcPhaser = vfcPhaser(:);
    vfcReceiveSignal = vfcTransmitSignal.*vfcPhaser;

    % sample offset
    vfcFFT = fft(vfcReceiveSignal);
    vOmega = 2*pi*( 0:length(vfcFFT)-1 ).'/length(vfcFFT);
    vfcReceiveSignal = ifft(vfcFFT.*fftshift(exp(1j*vOmega*stChannel.fSampleOffset)));

    % AWGN Channel
    vfcReceiveSignal = awgn(vfcReceiveSignal,stChannel.fSNRdB);

    %Multipath channel
    vfcReceiveSignal = conv(vfcReceiveSignal,stChannel.vfImpulseResponse);

    % Multipath channel
    if length(stChannel.vfImpulseResponse) > 1
        vfcReceiveSignal = conv(vfcReceiveSignal, stChannel.vfImpulseResponse, 'same');
    end
end
