function [iFrame,iBin] = frame_detect(Rlk,Plk,stOFDM)


lkMetrik = abs(xcorr2(Rlk,Plk));
%meshc(lkMetrik)

[value iFrame] = max(max(lkMetrik,[],2));
[value iBin] = max(max(lkMetrik,[],1));
iBin = iBin-stOFDM.iNfft;