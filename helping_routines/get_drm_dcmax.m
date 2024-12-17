%**************************************************************************
% Copyright:     (c) Univerity of Applied Science Rosenheim
% MODULE:         maximum used carrier index below DC-carrier
% PROJECT:       DRM Project
% ABBREVIATION:	 drm
% COMPILER:      MATLAB 2006a
% LANGUAGE:      Matlab Interpreter.
% AUTHOR:        Prof. Dr. Markus Stichler
% HISTORY:		 Initial version 11.05.2007
%**************************************************************************
%   returns maximum used carrier index below DC-carrier
%   for DRM Robustness Mode MODE
%   References:
%     [1] ETSI: Digital Radio Mondiale System Specification
%**************************************************************************
% function call:
%------------------------
% dcmax = get_drm_dcmax(mode);
%
% Output:
%------------------------
% dcmax: maximum used carrier index below DC-carrier
%
% Input:
%------------------------
% MODE       : DRM Robustness Mode (A = 1, B = 2, C = 3, D = 4)
%**************************************************************************
function dcmax = get_drm_dcmax(mode)

switch(mode)
    case 1
        dcmax = -2;
    case 2
        dcmax = -1;
    case 3
        dcmax = -1;
    case 4
        dcmax = -1;
    otherwise
        error('invalid DRM robustness mode !');
end

