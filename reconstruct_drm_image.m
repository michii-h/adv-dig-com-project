function [reconstructed_image, received_call_sign] = reconstruct_drm_image(Rlk, stDRM, M, image_size)
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

    % Template for the extraction of the data
    dataTemplate = repmat(get_drm_data_template_frame(stDRM.mode, stDRM.occupancy), size(Rlk,1)/15, 1);

    % Initialize empty vector for the data
    receivedSymbols = [];

    % Fill the vector row by row with the data from Rlk
    for iRow = 1:size(Rlk, 1)
        receivedSymbols = [receivedSymbols, Rlk(iRow, dataTemplate(iRow, :) == 1)];
    end

    % QAM-Demodulation
    demodulatedData = qamdemod(receivedSymbols', M);

    % Convert data to binary (4 bits per symbol)
    binaryData4 = de2bi(demodulatedData, 4, 'left-msb');

    % Reshape data to 8-bit per row
    [m, n] = size(binaryData4);
    binaryData8 = reshape(binaryData4', n*2, m/2)';

    % Reconversion to integer
    receivedData = bi2de(binaryData8, 'left-msb');

    % Extract call sign (first 6 bytes)
    received_call_sign = receivedData(1:6);

    % Extract image data
    imageData = receivedData(7:prod(image_size)+6);

    % Reconstruct image from received data
    reconstructed_image = reshape(imageData, image_size);
end
