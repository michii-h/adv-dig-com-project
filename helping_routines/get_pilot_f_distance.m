%**************************************************************************
% Copyright:     (c) Univerity of Applied Science Rosenheim
% MODULE:        gain pilot distance on carrier axis
% PROJECT:       DRM Project
% ABBREVIATION:	 drm
% COMPILER:      MATLAB 2006a
% LANGUAGE:      Matlab Interpreter.
% AUTHOR:        Prof. Dr. Markus Stichler
% HISTORY:		 Initial version 11.05.2007
%**************************************************************************
% returns gain pilot distance on carrier axis for DRM Robustness Mode MODE
%   References:
%     [1] ETSI: Digital Radio Mondiale System Specification
%**************************************************************************
% function call:
%------------------------
% pilot_f_distance = get_pilot_f_distance(mode);
%
% Output:
%------------------------
% pilot_f_distance: gain pilot distance on carrier axis
%
% Input:
%------------------------
% MODE       : DRM Robustness Mode (A = 1, B = 2, C = 3, D = 4)
%**************************************************************************


function pilot_f_distance = get_pilot_f_distance(mode)


switch(mode)
    case 1
        pilot_f_distance = 4;
    case 2
        pilot_f_distance = 2;
    case 3
        pilot_f_distance = 2;
    case 4
        pilot_f_distance = 1;
    otherwise
        error('Invalid DRM Robustness Mode !');
end

