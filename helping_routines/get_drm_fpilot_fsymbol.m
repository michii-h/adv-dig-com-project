%**************************************************************************
% Copyright:     (c) Univerity of Applied Science Rosenheim
% MODULE:        frequency pilot symbol in frequency domain
% PROJECT:       DRM Project
% ABBREVIATION:	 drm
% COMPILER:      MATLAB 2006a
% LANGUAGE:      Matlab Interpreter.
% AUTHOR:        Prof. Dr. Markus Stichler
% HISTORY:		 Initial version 11.05.2007
%**************************************************************************
%   returns frequency pilot symbol in frequency domain for DRM Robustness Mode MODE and spectrum occupancy
%   OCCUPANCY
%   References:
%     [1] ETSI: Digital Radio Mondiale System Specification
%**************************************************************************
% function call:
%------------------------
% fpilot_fsymbol = get_drm_fpilot_fsymbol(mode,occupancy);
%
% Output:
%------------------------
% fpilot_fsymbol: frequency pilot symbol in frequency domain
%
% Input:
%------------------------
% MODE       : DRM Robustness Mode (A = 1, B = 2, C = 3, D = 4)
% OCCUPANCY  : spectrum occupancy 0 - 5
%**************************************************************************
function fpilot_fsymbol = get_drm_fpilot_fsymbol(mode,occupancy)

fpilot_position = get_drm_fpilot_position(mode);
fpilot_data = get_drm_fpilot_data(mode);
fpilot_fsymbol = zeros(1,get_drm_n_useful(mode,occupancy));
dc_position = get_drm_dc_position(mode,occupancy);

for i = 1:length(fpilot_position)
    fpilot_fsymbol(fpilot_position(i) + dc_position) = fpilot_data(i);
end
