function [reconstructed_image, received_call_sign] = reconstruct_drm_image(Rlk, stDRM, M, image_size, expected_callsign)
% RECONSTRUCT_DRM_IMAGE Extracts and reconstructs image data from received DRM frame
%
% Inputs:
%   Rlk        - Received and equalized DRM frame
%   stDRM      - Structure containing DRM parameters (mode and occupancy)
%   M          - QAM modulation order (typically 16)
%   image_size - Size of the original image [rows, columns, channels]
%
% Outputs:
%   reconstructed_image - The reconstructed image
%   received_call_sign  - The extracted call sign characters
%
% Example:
%   [image, call_sign] = reconstruct_drm_image(Rlk, stDRM, 16, [480 640 3]);

    % Number of Symbols per frame
    iNOfSymbols = get_drm_symbols_per_frame(stDRM.mode);

    % Template for the extraction of the data
    dataTemplate = repmat(get_drm_data_template_frame(stDRM.mode, stDRM.occupancy), size(Rlk,1)/iNOfSymbols, 1);

    % Initialize empty vector for the data
    receivedSymbols = [];

    % Fill the vector row by row with the data from Rlk
    for iRow = 1:size(Rlk, 1)
        receivedSymbols = [receivedSymbols, Rlk(iRow, dataTemplate(iRow, :) == 1)];
    end

    % QAM-Demodulation
    %demodulatedData = qamdemod(receivedSymbols', M);
    if M == 4
        demodulatedData = qamdemod(receivedSymbols',M,'UnitAveragePower',true);
    else
        demodulatedData = qamdemod(receivedSymbols',M);
    end

    % Convert data to binary (4 bits per symbol)
    binaryData4 = de2bi(demodulatedData, sqrt(M), 'left-msb');

    % Cut off padding (for full frame) at end
    iNOfBits = (prod(image_size)+6)*8;
    iLenData = iNOfBits/sqrt(M);
    binaryData4 = binaryData4(1:iLenData,:);

    rescale_factor = 8/sqrt(M);

    % Reshape data to 8-bit per row
    [m, n] = size(binaryData4);
    %binaryData8 = reshape(binaryData4', n*rescale_factor, m/rescale_factor)';
    binaryData8 = reshape(binaryData4', n*rescale_factor, [])';

    % Reconversion to integer
    receivedData = bi2de(binaryData8, 'left-msb');

    % Find call sign
    start_index_callsign = strfind(receivedData', uint8(double(expected_callsign)));

    % Reshape Data
    if size(start_index_callsign) ~= 0
        % receivedData = [receivedData(start_index_callsign:end); receivedData(1:start_index_callsign)];
        receivedData = circshift(receivedData, start_index_callsign-1);
    else
        warning('Callsign not found');
    end

    % Extract call sign (first 6 bytes)
    received_call_sign = receivedData(1:6);

    % Extract image data
    imageData = receivedData(7:prod(image_size)+6);

    % Reconstruct image from received data
    reconstructed_image = reshape(imageData, image_size);
end
