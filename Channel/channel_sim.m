function vfcReceiveSignal = channel_sim(vfcTransmitSignal,stChannel)
    % Simulate iSampleShift = 17 between Tx and Rx
    iSampleShift = 0;
    vfcTransmitSignal = circshift(vfcTransmitSignal,[iSampleShift 0]);


    % Frequency- and Phase-Offset
    vfcPhaser = exp(j*stChannel.fOmegaOffset*[0:length(vfcTransmitSignal)-1]+j*stChannel.fPhaseOffset);
    vfcPhaser = vfcPhaser(:);
    vfcReceiveSignal = vfcTransmitSignal.*vfcPhaser;

    % sample offset
    vfcFFT = fft(vfcReceiveSignal);
    vOmega = 2*pi*[0:length(vfcFFT)-1].'/length(vfcFFT);
    vfcReceiveSignal = ifft(vfcFFT.*fftshift(exp(j*vOmega*stChannel.fSampleOffset)));

    % AWGN Channel
    vfcReceiveSignal = awgn(vfcReceiveSignal,stChannel.fSNRdB);

    %Multipath channel
    vfcReceiveSignal = conv(vfcReceiveSignal,stChannel.vfImpulseResponse);
end
