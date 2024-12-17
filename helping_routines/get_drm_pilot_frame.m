%**************************************************************************
% Copyright:     (c) Univerity of Applied Science Rosenheim
% MODULE:        complete drm frame
% PROJECT:       DRM Project
% ABBREVIATION:	 drm
% COMPILER:      MATLAB 2006a
% LANGUAGE:      Matlab Interpreter.
% AUTHOR:        Prof. Dr. Markus Stichler
% HISTORY:		 Initial version 11.05.2007
%**************************************************************************
%   returns DRM pilot frame for DRM Robustness Mode MODE and spectrum occupancy OCCUPANCY
%   References:
%     [1] ETSI: Digital Radio Mondiale System Specification
%**************************************************************************
% function call:
%------------------------
% pilot_frame = get_drm_pilot_frame(mode,occupancy);
%
% Output:
%------------------------
% pilot_frame: frame with complex pilot data vector
%
% Input:
%------------------------
% MODE       : DRM Robustness Mode (A = 1, B = 2, C = 3, D = 4)
% OCCUPANCY  : spectrum occupancy 0 - 5
%**************************************************************************
function pilot_frame = get_drm_pilot_frame(mode,occupancy)
% Initailization
% load pilots
fpilot_fsymbol = get_drm_fpilot_fsymbol(mode,occupancy);
tpilot_fsymbol = get_drm_tpilot_fsymbol(mode,occupancy);
[gpilot_position,gpilot_data] = get_drm_gpilot_position_data(mode,occupancy);
boosted_position =  get_drm_boosted_position(mode,occupancy);

symbols_per_frame = get_drm_symbols_per_frame(mode);
n_useful = get_drm_n_useful(mode,occupancy);
dc_position = get_drm_dc_position(mode,occupancy);

boost_factor = sqrt(2);

pilot_frame = zeros(symbols_per_frame, n_useful);

% insert frequency pilots
for l = 1:symbols_per_frame
    pilot_frame(l,:) = fpilot_fsymbol;
end

% special case Mode D
if (mode == 4)
    special_fpilot = get_drm_special_fpilot;
    for l = 1:symbols_per_frame/2
        for i = 1:length(special_fpilot)
            pilot_frame(2*l-1,special_fpilot(i) + dc_position) = -pilot_frame(2*l-1,special_fpilot(i) + dc_position);
        end
    end
end

% insert timing pilots
pilot_frame(1,:) = pilot_frame(1,:) + tpilot_fsymbol;

% insert gain pilots
for l = 1:symbols_per_frame
    for i = 1:length(gpilot_position(l,:))
        if (pilot_frame(l,gpilot_position(l,i) + dc_position) == 0)
            pilot_frame(l,gpilot_position(l,i) + dc_position) = gpilot_data(l,i);
        end
    end
    % boost outer pilots
    for i = 1:length(boosted_position)
        pilot_frame(l,boosted_position(i) + dc_position) = boost_factor * pilot_frame(l,boosted_position(i) + dc_position);
    end
end