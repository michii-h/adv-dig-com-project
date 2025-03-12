function [txPower, margin, snr] = calculate_link_budget(stSat)
  % Common calculations
  lambda = stSat.c / stSat.uplinkFreq;

  % Define parameters for both calculation methods
  antennaEfficiency = 0.5;      % 50% efficiency as used in example
  cableLength = 12;             % meters of cable
  cableLoss = 0.3;              % dB/m for Ecoflex 7mm

  % Antenna gain calculation
  antennaGain = 10*log10(((pi*stSat.antennaDiam)/lambda)^2 * antennaEfficiency);

  % Free space path loss
  fsl = 20*log10(4*pi*stSat.slantRange/lambda);

  % Cable loss calculation
  totalCableLoss = cableLength * cableLoss;

  % System noise calculation
  noiseTemp = 3;
  systemTemp = noiseTemp / 10^(stSat.gt/10);
  noiseFloor = 10*log10(stSat.k * systemTemp * stSat.bandwidth);

  % noiseFloor = 10*log10(stSat.k * stSat.transponderBW) + stSat.gt;

  % Calculate required received power at satellite
  requiredRxPower = noiseFloor + stSat.targetSNR;

  % Calculate required transmit power (dBW)
  requiredTxPowerDB = requiredRxPower + fsl + totalCableLoss - antennaGain;

  % Convert to linear power (Watts)
  txPower = 10^(requiredTxPowerDB/10);

  % Calculate actual SNR with this power
  rxPowerAtSat = 10*log10(txPower) + antennaGain - fsl - totalCableLoss;
  snr = rxPowerAtSat - noiseFloor;
  margin = snr - stSat.targetSNR;

  % ---- Print results ----
  fprintf('Link budget calculations:\n');
  fprintf(' Free space path loss: %.2f dB\n', fsl);
  fprintf(' Antenna gain: %.2f dBi\n', antennaGain);
  fprintf(' Cable loss: %.2f dB\n', totalCableLoss);
  fprintf(' Noise floor: %.2f dBW\n', noiseFloor);
end
