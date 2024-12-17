%**************************************************************************
% Copyright:     (c) Univerity of Applied Science Rosenheim
% MODULE:        frequency pilot data for DRM Robustness Mode MODE
% PROJECT:       DRM Project
% ABBREVIATION:	 drm
% COMPILER:      MATLAB 2006a
% LANGUAGE:      Matlab Interpreter.
% AUTHOR:        Prof. Dr. Markus Stichler
% HISTORY:		 Initial version 11.05.2007
%**************************************************************************
%   returns frequency pilot data for DRM Robustness Mode MODE
%   References:
%     [1] ETSI: Digital Radio Mondiale System Specification
%**************************************************************************
% function call:
%------------------------
% fpilot_data = get_drm_fpilot_data(mode);
%
% Output:
%------------------------
% fpilot_data: pilot data for DRM Robustness Mode MODE
%
% Input:
%------------------------
% MODE       : DRM Robustness Mode (A = 1, B = 2, C = 3, D = 4)
%**************************************************************************
function fpilot_data = get_drm_fpilot_data(mode)

magnitude = sqrt(2);
increment = pi / 512;

switch(mode)
    case 1
        fpilot_data = [
            magnitude * exp(j * 205 * increment)
            magnitude * exp(j * 836 * increment)
            magnitude * exp(j * 215 * increment)
            ];
    case 2
        fpilot_data = [
            magnitude * exp(j * 331 * increment)
            magnitude * exp(j * 651 * increment)
            magnitude * exp(j * 555 * increment)
            ];
    case 3
        fpilot_data = [
            magnitude * exp(j * 214 * increment)
            magnitude * exp(j * 392 * increment)
            magnitude * exp(j * 242 * increment)
            ];
    case 4
        fpilot_data = [
            magnitude * exp(j * 788 * increment)
            magnitude * exp(j *1014 * increment)
            magnitude * exp(j * 332 * increment)
            ];
    otherwise
        error('invalid DRM robustness mode !');
end


