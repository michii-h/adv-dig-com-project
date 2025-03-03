function [Rlk_corrected, freq_offset, phase_offset] = optimize_fine_sync(Rlk, Plk, iNfft)
% OPTIMIZE_FINE_SYNC Optimizes fine synchronization using FFT implementation
%   This function performs fine synchronization for frequency and phase correction
%   using an FFT-based algorithm to increase SNR performance
%
%   Parameters:
%   Rlk - Received OFDM symbols after FFT
%   Plk - Pilot pattern for the frame
%   iNfft - FFT size
%
%   Returns:
%   Rlk_corrected - Frequency/phase corrected OFDM symbols
%   freq_offset - Estimated frequency offset
%   phase_offset - Estimated phase offset

% Initialize outputs
Rlk_corrected = Rlk;
freq_offset = 0;
phase_offset = 0;

% Get pilot positions
[rows, cols] = find(Plk ~= 0);
pilot_positions = [rows cols];

% Extract pilot values
pilot_values = zeros(size(pilot_positions, 1), 1);
for i = 1:size(pilot_positions, 1)
    r = pilot_positions(i, 1);
    c = pilot_positions(i, 2);
    pilot_values(i) = Rlk(r, c) / Plk(r, c);
end

% 1. Frequency offset estimation using FFT-based algorithm
% Group pilots by symbol (row)
unique_rows = unique(rows);
num_symbols = length(unique_rows);

% For each symbol, compute frequency offset
freq_offsets = zeros(num_symbols, 1);
for i = 1:num_symbols
    symbol_idx = rows == unique_rows(i);
    symbol_pilots = pilot_positions(symbol_idx, :);
    symbol_values = pilot_values(symbol_idx);

    % Extract pilot positions and values for this symbol
    k = symbol_pilots(:, 2) - iNfft/2 - 1; % Centered subcarrier index
    phases = unwrap(angle(symbol_values));

    % Use FFT-based algorithm for better frequency estimation
    max_offset = 0.5; % Maximum frequency offset (normalized)
    fft_size = 1024;
    correlation = zeros(fft_size, 1);

    % Create grid for frequency offsets
    freq_grid = linspace(-max_offset, max_offset, fft_size);

    % Calculate correlation for each potential frequency offset
    for f = 1:fft_size
        test_freq = freq_grid(f);
        exp_values = exp(-1j * 2 * pi * test_freq * k);
        correlation(f) = abs(sum(exp_values .* symbol_values));
    end

    % Find the frequency offset with maximum correlation
    [~, max_idx] = max(correlation);
    freq_offsets(i) = freq_grid(max_idx);
end

% Average frequency offset across all symbols
freq_offset = mean(freq_offsets);

% 2. Apply frequency correction
k_indices = (0:size(Rlk, 2)-1) - iNfft/2 - 1;
for l = 1:size(Rlk, 1)
    phase_correction = exp(-1j * 2 * pi * freq_offset * k_indices);
    Rlk_corrected(l, :) = Rlk(l, :) .* phase_correction;
end

% 3. Phase offset estimation after frequency correction
pilot_values_corrected = zeros(size(pilot_positions, 1), 1);
for i = 1:size(pilot_positions, 1)
    r = pilot_positions(i, 1);
    c = pilot_positions(i, 2);
    pilot_values_corrected(i) = Rlk_corrected(r, c) / Plk(r, c);
end

phase_offset = angle(mean(pilot_values_corrected));

% 4. Apply phase correction
phase_correction = exp(-1j * phase_offset);
Rlk_corrected = Rlk_corrected * phase_correction;

% 5. Residual frequency offset estimation per symbol
for l = 1:size(Rlk_corrected, 1)
    symbol_idx = rows == l;
    if sum(symbol_idx) < 2
        continue; % Need at least 2 pilots
    end

    symbol_pilots = pilot_positions(symbol_idx, :);
    symbol_values = zeros(sum(symbol_idx), 1);

    for i = 1:length(symbol_values)
        r = symbol_pilots(i, 1);
        c = symbol_pilots(i, 2);
        symbol_values(i) = Rlk_corrected(r, c) / Plk(r, c);
    end

    k = symbol_pilots(:, 2) - iNfft/2 - 1;
    phases = unwrap(angle(symbol_values));

    % Linear regression for phase vs. subcarrier index
    if length(k) > 1
        p = polyfit(k, phases, 1);
        residual_freq = p(1) / (2 * pi);

        % Apply residual frequency correction for this symbol
        residual_correction = exp(-1j * 2 * pi * residual_freq * k_indices);
        Rlk_corrected(l, :) = Rlk_corrected(l, :) .* residual_correction;
    end
end

% Print results
fprintf('Fine Synchronization Results:\n');
fprintf('  Frequency Offset: %.6f (normalized)\n', freq_offset);
fprintf('  Phase Offset: %.6f radians\n', phase_offset);

end
