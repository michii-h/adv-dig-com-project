%**************************************************************************
% Copyright:     (c) Univerity of Applied Science Rosenheim
% MODULE:        number of symbols per frame
% PROJECT:       DRM Project
% ABBREVIATION:	 drm
% COMPILER:      MATLAB 2006a
% LANGUAGE:      Matlab Interpreter.
% AUTHOR:        Prof. Dr. Markus Stichler
% HISTORY:		 Initial version 11.05.2007
%**************************************************************************
%   returns number of symbols per frame for DRM Robustness Mode MODE
%   References:
%     [1] ETSI: Digital Radio Mondiale System Specification
%**************************************************************************
% function call:
%------------------------
% symbols_per_frame = get_drm_symbols_per_frame(mode);
%
% Output:
%------------------------
% symbols_per_frame: number of symbols per frame
%
% Input:
%------------------------
% MODE       : DRM Robustness Mode (A = 1, B = 2, C = 3, D = 4)
%**************************************************************************


function symbols_per_frame = get_drm_symbols_per_frame(mode)

switch(mode)
    case 1
        symbols_per_frame = 15;
    case 2
        symbols_per_frame = 15;
    case 3
        symbols_per_frame = 20;
    case 4
        symbols_per_frame = 24;
    otherwise
        error('invalid DRM robustness mode !');
end
