%**************************************************************************
% Copyright:     (c) Univerity of Applied Science Rosenheim
% MODULE:        timing pilot symbol in frequency domain
% PROJECT:       DRM Project
% ABBREVIATION:	 drm
% COMPILER:      MATLAB 2006a
% LANGUAGE:      Matlab Interpreter.
% AUTHOR:        Prof. Dr. Markus Stichler
% HISTORY:		 Initial version 11.05.2007
%**************************************************************************
%   returns timing pilot symbol in frequency domain for DRM Robustness Mode MODE and spectrum occupancy
%   OCCUPANCY
%   References:
%     [1] ETSI: Digital Radio Mondiale System Specification
%**************************************************************************
% function call:
%------------------------
% tpilot_fsymbol = get_drm_tpilot_fsymbol(mode,occupancy);
%
% Output:
%------------------------
% tpilot_fsymbol: timing pilot symbol in frequency domain
%
% Input:
%------------------------
% MODE       : DRM Robustness Mode (A = 1, B = 2, C = 3, D = 4)
% OCCUPANCY  : spectrum occupancy 0 - 5
%**************************************************************************

function tpilot_fsymbol = get_drm_tpilot_fsymbol(mode,occupancy)

tpilot_position = get_drm_tpilot_position(mode);
tpilot_data = get_drm_tpilot_data(mode);
tpilot_fsymbol = zeros(1,get_drm_n_useful(mode,occupancy));
dc_position = get_drm_dc_position(mode,occupancy);

for i = 1:length(tpilot_position)
    tpilot_fsymbol(tpilot_position(i) + dc_position) = tpilot_data(i);
end

