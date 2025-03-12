function check_drm_bandwidth(fs, mode, occ)
    fprintf('DRM Bandwidth for 2.7kHz constraint:\n');

    % Get DRM parameters
    iNfft = get_drm_n_useful(mode, occ);
    kmin = get_drm_kmin(mode, occ);
    kmax = get_drm_kmax(mode, occ);

    % Calculate bandwidth
    subcarrierSpacing = fs / iNfft;
    carriers = kmax - kmin + 1;
    bandwidth = carriers * subcarrierSpacing / 1000; % in kHz

    if (bandwidth <= 2.7); sign = '✓' ; else; sign = '✗'; end

    fprintf('  Mode: %d\n', mode);
    fprintf('  Occupancy: %d\n', occ);
    fprintf('  FFT Size: %d\n', iNfft);
    fprintf('  Carriers: %d\n', carriers);
    fprintf('  Spacing: %.2f Hz\n', subcarrierSpacing);
    fprintf('  Bandwidth: %.2f kHz %s\n', bandwidth, sign);
end
