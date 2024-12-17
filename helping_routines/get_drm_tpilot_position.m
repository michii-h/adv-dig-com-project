%**************************************************************************
% Copyright:     (c) Univerity of Applied Science Rosenheim
% MODULE:        timing pilot positions
% PROJECT:       DRM Project
% ABBREVIATION:	 drm
% COMPILER:      MATLAB 2006a
% LANGUAGE:      Matlab Interpreter.
% AUTHOR:        Prof. Dr. Markus Stichler
% HISTORY:		 Initial version 11.05.2007
%**************************************************************************
%  returns timing pilot positions for DRM Robustness Mode MODE
%   References:
%     [1] ETSI: Digital Radio Mondiale System Specification
%**************************************************************************
% function call:
%------------------------
% tpilot_position = get_drm_tpilot_position(mode);
%
% Output:
%------------------------
% tpilot_position: timing pilot positions
%
% Input:
%------------------------
% MODE       : DRM Robustness Mode (A = 1, B = 2, C = 3, D = 4)
%**************************************************************************


function tpilot_position = get_drm_tpilot_position(mode)


switch(mode)
    case 1
        tpilot_position = [17 19 21 28 29 32 33 39 40 41 53 55 56 60 61 63 71 73];
    case 2
        tpilot_position = [14 18 20 24 26 32 36 42 44 49 50 54 56 62 66 68];
    case 3
        tpilot_position = [ 8 10 12 14 16 18 22 24 28 30 32 36 38 42 45 46];
    case 4
        tpilot_position = [ 5  6  8  9 11 12 14 15 17 18 20 23 24 26 27 29 30 32];
    otherwise
        error('invalid DRM robustness mode !');
end
