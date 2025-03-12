function calculate_drm_bandwidth(fs, mode, occ)
    fprintf('DRM Bandwidth for 2.7kHz constraint:\n');

    fprintf('---------------------------------------------------\n');
    fprintf('Mode | Occ | FFT Size | Carriers | Spacing (Hz) | BW (kHz)\n');
    fprintf('---------------------------------------------------\n');

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
end
