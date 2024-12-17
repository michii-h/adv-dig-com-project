%**************************************************************************
% Copyright:     (c) Univerity of Applied Science Rosenheim
% MODULE:        frequency pilot positions
% PROJECT:       DRM Project
% ABBREVIATION:	 drm
% COMPILER:      MATLAB 2006a
% LANGUAGE:      Matlab Interpreter.
% AUTHOR:        Prof. Dr. Markus Stichler
% HISTORY:		 Initial version 11.05.2007
%**************************************************************************
%   returns frequency pilot positions for DRM Robustness Mode MODE
%   References:
%     [1] ETSI: Digital Radio Mondiale System Specification
%**************************************************************************
% function call:
%------------------------
% fpilot_position = get_drm_fpilot_position(mode);
%
% Output:
%------------------------
% fpilot_position: frequency pilot positions
%
% Input:
%------------------------
% MODE       : DRM Robustness Mode (A = 1, B = 2, C = 3, D = 4)
%**************************************************************************
function fpilot_position = get_drm_fpilot_position(mode)

switch(mode)
    case 1
        fpilot_position = [18 54 72];
    case 2
        fpilot_position = [16 48 64];
    case 3
        fpilot_position = [11 33 44];
    case 4
        fpilot_position = [7 21 28];
    otherwise
        error('invalid DRM robustness mode !');
end
