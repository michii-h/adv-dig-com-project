%**************************************************************************
% Copyright:     (c) Univerity of Applied Science Rosenheim
% MODULE:        minimum used carrier index below DC-carrier
% PROJECT:       DRM Project
% ABBREVIATION:	 drm
% COMPILER:      MATLAB 2006a
% LANGUAGE:      Matlab Interpreter.
% AUTHOR:        Prof. Dr. Markus Stichler
% HISTORY:		 Initial version 11.05.2007
%**************************************************************************
%   returns minimum used carrier index below DC-carrier
%   for DRM Robustness Mode MODE
%   References:
%     [1] ETSI: Digital Radio Mondiale System Specification
%**************************************************************************
% function call:
%------------------------
% dcmin = get_drm_dcmin(mode);
%
% Output:
%------------------------
% dcmin: minimum used carrier index below DC-carrier
%
% Input:
%------------------------
% MODE       : DRM Robustness Mode (A = 1, B = 2, C = 3, D = 4)
%**************************************************************************

function dcmin = get_drm_dcmin(mode)

switch(mode)
    case 1
        dcmin = 2;
    case 2
        dcmin = 1;
    case 3
        dcmin = 1;
    case 4
        dcmin = 1;
    otherwise
        error('invalid DRM robustness mode !');
end
