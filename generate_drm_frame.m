function [Slk, M, image_size, iNofFramesNeeded, iNOfFrames] = generate_drm_frame(stDRM, image_path, call_sign_str)
% GENERATE_DRM_FRAME Generates a DRM frame with embedded image and call sign
%   [Slk, M, image_size, iNofFramesNeeded, iNOfFrames] = generate_drm_frame(stDRM, image_path, call_sign_str)
%
%   Inputs:
%     stDRM - Structure containing DRM parameters (mode, occupancy)
%     image_path - Path to the image file to be embedded
%     call_sign_str - String containing the call sign
%
%   Outputs:
%     Slk - The generated DRM frame
%     M - Modulation order
%     image_size - Size of the original image for reconstruction
%     iNofFramesNeeded - Number of frames needed for the data
%     iNOfFrames - Number of frames

% Number of Symbols per frame
iNOfSymbols = get_drm_symbols_per_frame(stDRM.mode);

% Initialize Frame Matrix
Slk = get_drm_data_template_frame(stDRM.mode, stDRM.occupancy);

% Load Image
image = imread(image_path);
image_size = size(image);

% Reconstruct the image to a vector
viImage = image(:)';

% Concat Call Sign
call_sign = uint8(double(call_sign_str));
viData = [call_sign, viImage];

%convert to binary 8-bits per integer
binaryData8 = de2bi(viData, 8, 'left-msb');

% reshape to 4-bits per row
[m,n]=size(binaryData8);
binaryData4 = reshape(binaryData8', n/2, m*2)';

% converte back to integer
viData = bi2de(binaryData4,'left-msb');

% Datamapping: M-QAM
M = 16;     % 4 bits per row -> 16 possible values
% viDlk = qammod(viData,M, 'UnitAveragePower',true);
viDlk = qammod(viData,M);

% concat rows one after another
viDlk = reshape(viDlk.',1,[]);

% Set Data in Slk
vSlk = reshape(Slk',1,[]);

% Calculate Data per Frame
DataPerFrame = sum(vSlk);

% Calculate needed number of Frames
slkCtr = 1;
dataCtr = 1;
resultCtr = 1;
iNofFramesNeeded = ceil( size(viDlk, 2) / DataPerFrame );

% Preallocate viDataPadded
viDataPadded = zeros(1, iNofFramesNeeded * DataPerFrame);

% Pad Data with Zeros at Pilot Positions
while dataCtr <= size(viDlk, 2)
    if vSlk(slkCtr) == 1
        viDataPadded(resultCtr) = viDlk(dataCtr);
        dataCtr = dataCtr + 1;
    else
        % already zero
    end
    slkCtr = mod(slkCtr, 15*256) + 1; % 15*256 = len of Slk
    resultCtr = resultCtr + 1;
end

% Zeropadding at end for full DRM frames
nRows = ceil(size(viDataPadded, 2) / 256);
nRowsPad = 15 - mod(nRows, 15);
nRows = nRows + nRowsPad;

zData = zeros(1, nRows * 256);
zData(1:size(viDataPadded,2)) = viDataPadded;

% Reshape to Frameshape
viDataPadded = reshape(zData, 256, [])';
Slk = viDataPadded;

% Generate Pilots
Plk = get_drm_pilot_frame(stDRM.mode,stDRM.occupancy);

% Duplicate Pilots to length of Slk
iNOfFrames = nRows / iNOfSymbols;
Plk = repmat(Plk,[iNOfFrames 1]);

% Set Pilots in Slk
Slk(Plk ~= 0) = Plk(Plk ~= 0);

end
