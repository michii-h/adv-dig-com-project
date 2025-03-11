function vfcReceiveSignal = sync(vfcReceiveSignal, iNs, iNg, iNfft, SwitchDemoSync)
% Number of Symbols in CaptureBuffer
iNOfSymbolsTotal = floor(length(vfcReceiveSignal)/iNs);

% Timing Recovery with Cyclic Prefix
vIndexGI = [0:iNg-1].';
R = zeros(1,iNs);

for l = 0:iNOfSymbolsTotal-2
    for k = 1:iNs
        R(k) = R(k) + ...
            sum(conj(vfcReceiveSignal(k+vIndexGI+l*iNs)).*vfcReceiveSignal(k+vIndexGI+iNfft+l*iNs));
    end
end

[~ ,iStartSample] = max(abs(R));
dOmegaEst = angle(R(iStartSample))/iNfft;

% Compensate Frequency Offset
vfcPhaser = exp(-1j*dOmegaEst*[0:length(vfcReceiveSignal)-1]);
vfcPhaser = vfcPhaser(:);

% Adjust to Start Sample
vfcReceiveSignal = vfcReceiveSignal(iStartSample:end);
iNOfSymbolsTotal = floor(length(vfcReceiveSignal)/iNs);
vfcReceiveSignal = vfcReceiveSignal(1:iNOfSymbolsTotal*iNs);

if SwitchDemoSync
    figure(101)
    plot([0:iNs-1],abs(R))
    xlabel('k')
    ylabel('R_{xy}(k)')
    hold on
    plot(iStartSample, abs(R(iStartSample)),'ro')
    hold off

    % figure;
    % plot(abs(sum(reshape(filter(ones(1,iNg),1,vfcReceiveSignal.*conj(circshift(vfcReceiveSignal,[iNfft 0]))),iNs,iNOfSymbols),2)))
end
end