% Cyclic prefix metric to estimate
% a) Startsample of OFDM symbol: iStart in samples
% b) fractional frequency offset: dOmega in rad
function [dOmega, iStart] = timing_metric(vfcReceiveSignal,stOFDM,SwitchDemoSync)
%% Variante 1: Formel
% Initialization
cpMetric = zeros(1,stOFDM.iNs);
iNOfSymbols = floor(length(vfcReceiveSignal)/stOFDM.iNs);
iCPSamples = [0:stOFDM.iNg-1];
%tic
for l = 0:iNOfSymbols-2
    for n = 1:stOFDM.iNs
        
        cpMetric(n)=cpMetric(n) + ...
            vfcReceiveSignal(n+l*stOFDM.iNs+iCPSamples)'*vfcReceiveSignal(n+l*stOFDM.iNs+iCPSamples+stOFDM.iNfft) / ...
            (sqrt(vfcReceiveSignal(n+l*stOFDM.iNs+iCPSamples)'*vfcReceiveSignal(n+l*stOFDM.iNs+iCPSamples)) * ...
            sqrt(vfcReceiveSignal(n+l*stOFDM.iNs+iCPSamples+stOFDM.iNfft)'*vfcReceiveSignal(n+l*stOFDM.iNs+iCPSamples+stOFDM.iNfft)));
        
    end
end
%toc

cpMetric = cpMetric/iNOfSymbols;

[value iStart] = max(abs(cpMetric));
dOmega = angle(cpMetric(iStart))/stOFDM.iNfft;

if SwitchDemoSync
figure(10)
plot(abs(cpMetric))
grid
title('schmidl cox timing metric normalized as reliabilty factor')
xlabel('n \rightarrow')
ylabel('R(n) \rightarrow')
hold on
plot(iStart,value,'ro')
xline(iStart,'-',['start sample: n=' num2str(iStart)  ' / Reliabilty ' num2str(abs(cpMetric(iStart))*100) ' %']);
hold off
end


% %% Variante 2: Schnelle implementierung
% iNOfSymbols = floor(length(vfcReceiveSignal)/stOFDM.iNs);
% tic
% metric = (vfcReceiveSignal) .* conj(circshift(vfcReceiveSignal,[stOFDM.iNfft 0]));
% metric = filter(ones(stOFDM.iNg,1),1,metric);
% metric = metric(1:iNOfSymbols*stOFDM.iNs);
% mmetric = reshape(metric,stOFDM.iNs,iNOfSymbols);
% cpmetric = mean(mmetric,2);
% toc
% [value iStart] = max(abs(cpmetric));
% dOmega = angle(cpmetric(iStart));
