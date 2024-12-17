%**************************************************************************
% Copyright:     (c) Univerity of Applied Science Rosenheim
% MODULE:        dc position
% PROJECT:       DRM Project
% ABBREVIATION:	 drm
% COMPILER:      MATLAB 2006a
% LANGUAGE:      Matlab Interpreter.
% AUTHOR:        Prof. Dr. Markus Stichler
% HISTORY:		 Initial version 11.05.2007
%**************************************************************************
%   returns position of DC-carrier for DRM Robustness Mode MODE and
%   spectrum occupancy OCCUPANCY
%   References:
%     [1] ETSI: Digital Radio Mondiale System Specification
%**************************************************************************
% function call:
%------------------------
% dc_position = get_drm_dc_position(mode,occupancy)
%
% Output:
%------------------------
% dc_position
%
% Input:
%------------------------
% MODE       : DRM Robustness Mode (A = 1, B = 2, C = 3, D = 4)
% OCCUPANCY  : spectrum occupancy 0 - 5
%**************************************************************************
function dc_position = get_drm_dc_position(mode,occupancy)

switch(occupancy)
    case {0,1}
        dc_position = 1;
    case {2,3}
        dc_position = get_drm_n_useful(mode,occupancy) / 2 + 1;
    case {4,5}
        dc_position = get_drm_n_useful(mode,occupancy) / 4 + 1;
end

