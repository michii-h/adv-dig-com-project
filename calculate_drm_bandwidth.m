function calculate_drm_bandwidth(fs)
    % Define modes to check
    modes = [1, 2, 3, 4]; % Modes A, B, C, D
    occupancies = [0, 1, 2, 3, 4, 5]; % Different occupancies

    fprintf('DRM Bandwidth Calculations for 2.7kHz constraint:\n');
    fprintf('---------------------------------------------------\n');
    fprintf('Mode | Occ | FFT Size | Carriers | Spacing (Hz) | BW (kHz)\n');
    fprintf('---------------------------------------------------\n');

    for mode = modes
        for occ = occupancies
            try
                % Get DRM parameters
                iNfft = get_drm_n_useful(mode, occ);
                kmin = get_drm_kmin(mode, occ);
                kmax = get_drm_kmax(mode, occ);

                % Calculate bandwidth
                subcarrierSpacing = fs / iNfft;
                carriers = kmax - kmin + 1;
                bandwidth = carriers * subcarrierSpacing / 1000; % in kHz

                if (bandwidth <= 2.7); sign = '✓' ; else; sign = '✗'; end

                % Print results
                fprintf(' %d   | %d   | %4d     | %4d     | %7.2f    | %5.2f %s\n', ...
                    mode, occ, iNfft, carriers, subcarrierSpacing, bandwidth, sign);
            catch
                % Skip invalid combinations
            end
        end
    end
end
