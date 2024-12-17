%**************************************************************************
% Copyright:     (c) Univerity of Applied Science Rosenheim
% MODULE:        gain pilot frame with 1/x power in frequency domain
% PROJECT:       DRM Project
% ABBREVIATION:	 drm
% COMPILER:      MATLAB 2006a
% LANGUAGE:      Matlab Interpreter.
% AUTHOR:        Prof. Dr. Markus Stichler
% HISTORY:		 Initial version 11.05.2007
%**************************************************************************
%   returns gain pilot frame with 1/x power in frequency domain for DRM Robustness Mode MODE and
%   spectrum occupancy OCCUPANCY
%   References:
%     [1] ETSI: Digital Radio Mondiale System Specification
%**************************************************************************
% function call:
%------------------------
% gpilot_fframe = get_drm_gpilot_equalizer_fframe(mode,occupancy);
%
% Output:
%------------------------
% gpilot_fframe: gain pilot frame with 1/x power in frequency domain
%
% Input:
%------------------------
% MODE       : DRM Robustness Mode (A = 1, B = 2, C = 3, D = 4)
% OCCUPANCY  : spectrum occupancy 0 - 5
%**************************************************************************

function gpilot_fframe = get_drm_gpilot_equalizer_fframe(mode,occupancy)

% load pilots
[gpilot_position,gpilot_data] = get_drm_gpilot_position_data(mode,occupancy);
boosted_position =  get_drm_boosted_position(mode,occupancy);

symbols_per_frame = get_drm_symbols_per_frame(mode);
n_useful = get_drm_n_useful(mode,occupancy);
dc_position = get_drm_dc_position(mode,occupancy);

boost_factor = sqrt(0.5);

gpilot_fframe = zeros(symbols_per_frame, n_useful);

% Insert gain pilots
for l = 1:symbols_per_frame
    for i = 1:length(gpilot_position(l,:))
        if (gpilot_fframe(l,gpilot_position(l,i) + dc_position) == 0)
            gpilot_fframe(l,gpilot_position(l,i) + dc_position) = gpilot_data(l,i) / 2;
        end
    end
    % boost outer pilots
    for i = 1:length(boosted_position)
        gpilot_fframe(l,boosted_position(i) + dc_position) = boost_factor * gpilot_fframe(l,boosted_position(i) + dc_position);
    end
end
