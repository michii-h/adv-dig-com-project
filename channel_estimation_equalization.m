function Rlk_out = channel_estimation_equalization(Rlk, Plk, iNfft, iNg, iModeEst, iOcc, method, SwitchDemoSync)
% CHANNEL_ESTIMATION_EQUALIZATION Performs channel estimation and equalization
%
% Inputs:
%   Rlk            - Received OFDM symbols in frequency domain
%   Plk            - Pilot matrix
%   iNfft          - FFT size
%   iNg            - Guard interval length
%   iModeEst       - Estimated DRM mode
%   iOcc           - Spectrum occupancy
%   method         - Interpolation method ('Spline' or 'Wiener')
%   SwitchDemoSync - Flag for additional plots
%
% Outputs:
%   Rlk_out - Equalized received symbols

% Default method if not specified
if nargin < 7
    method = 'Wiener';
end

% Default plot flag if not specified
if nargin < 8
    SwitchDemoSync = false;
end

% Get DRM parameters
dc = get_drm_dc_position(iModeEst, iOcc);
kmin = get_drm_kmin(iModeEst, iOcc) + dc;
kmax = get_drm_kmax(iModeEst, iOcc) + dc;
kInterpolate = kmin:kmax;

% Copy the input for output
Rlk_out = Rlk;

% Process each OFDM symbol
for l = 1:size(Plk, 1)
    % Get pilot positions for current symbol
    kIndex = find(Plk(l, :) ~= 0);

    % Get channel estimates at pilot positions
    Hk = Rlk(l, kIndex) ./ Plk(l, kIndex);

    % Initialize channel estimate vector
    Hint = zeros(1, iNfft);

    % Perform interpolation based on selected method
    switch method
        case 'Spline'
            Hint = spline_interpolation(kIndex, Hk, kInterpolate, iNfft);
        case 'Wiener'
            Hint = wiener_interpolation(kIndex, Hk, kInterpolate, iNfft, iNg, SwitchDemoSync);
        otherwise
            error('Unknown interpolation method: %s', method);
    end

    % Apply windowing to the channel estimate
    Hint(kInterpolate) = Hint(kInterpolate) .* hamming(1, length(kInterpolate))';

    % Calculate time domain impulse response (for analysis)
    h_est = ifft(fftshift(Hint));

    % Equalize the symbol
    Rlk_out(l, kInterpolate) = Rlk(l, kInterpolate) ./ Hint(kInterpolate);
end
end

function Hint = spline_interpolation(kIndex, Hk, kInterpolate, iNfft)
% SPLINE_INTERPOLATION Performs spline interpolation for channel estimation
%
% Simple cubic spline interpolation between pilot positions

Hint = zeros(1, iNfft);
Hint(kInterpolate) = interp1(kIndex, Hk, kInterpolate, 'spline');
end

function Hint = wiener_interpolation(kIndex, Hk, kInterpolate, iNfft, iNg, SwitchDemoSync)
% WIENER_INTERPOLATION Performs Wiener interpolation for channel estimation
%
% Wiener filtering for optimal MMSE channel estimation

% Initialize channel estimate vector
Hint = zeros(1, iNfft);

% Ensure Hk is a column vector
Hk = Hk(:);

% Calculate correlation matrix
mIndex = [];
for k = kIndex
    mIndex = [mIndex; kIndex - k];
end

% Set frequency domain correlation parameters
Omega_g = 0.5 * 1/iNg;  % Normalized Doppler frequency

% Calculate correlation matrix
R_ = sinc(mIndex * Omega_g);

% Add noise to correlation matrix for stability
SNR = 50;  % Assumed SNR in dB
N = 10^(-SNR/10);
R = R_ + eye(size(R_, 1)) * N;

% Pre-calculate inverse for efficiency
Rinv = inv(R);

% Perform interpolation for each subcarrier
for k = kInterpolate
    % Observation vector
    vIndex = kIndex - k;
    b = sinc(vIndex * Omega_g);
    b = b(:);

    % Wiener filter
    g = Rinv * b;

    % Apply filter to get channel estimate
    Hint(k) = g' * Hk;
end

% Display debug plots if enabled
if SwitchDemoSync
    figure(10);
    stem(abs(Hint));
    hold on;
    plot(kIndex, abs(Hk), 'ro');
    hold off;
    title('Channel Estimation');
    xlabel('Subcarrier Index');
    ylabel('Magnitude');
    legend('Interpolated', 'Pilots');
end
end
