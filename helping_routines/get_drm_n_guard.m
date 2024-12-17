%**************************************************************************
% Copyright:     (c) Univerity of Applied Science Rosenheim
% MODULE:        length of guard interval
% PROJECT:       DRM Project
% ABBREVIATION:	 drm
% COMPILER:      MATLAB 2006a
% LANGUAGE:      Matlab Interpreter.
% AUTHOR:        Prof. Dr. Markus Stichler
% HISTORY:		 Initial version 11.05.2007
%**************************************************************************
%   returns length of guard interval
%   for DRM Robustness Mode MODE and spectrum occupancy OCCUPANCY
%   References:
%     [1] ETSI: Digital Radio Mondiale System Specification
%**************************************************************************
% function call:
%------------------------
% n_guard = get_drm_n_guard(mode,occupancy);
%
% Output:
%------------------------
% n_guard: length of guard interval
%
% Input:
%------------------------
% MODE       : DRM Robustness Mode (A = 1, B = 2, C = 3, D = 4)
% OCCUPANCY  : spectrum occupancy 0 - 5
%**************************************************************************

function n_guard = get_drm_n_guard(mode,occupancy)


switch(mode)
    case 1
        switch(occupancy)
            case {0,1,2,3}
                n_guard = 32;
            case {4,5}
                n_guard = 64;
            otherwise
                error('invalid DRM spectrum occupancy !');
        end
    case 2
        switch(occupancy)
            case {0,1,2,3}
                n_guard = 64;
            case {4,5}
                n_guard = 128;
            otherwise
                error('invalid DRM spectrum occupancy !');
        end
    case 3
        switch(occupancy)
            case 3
                n_guard = 64;
            case 5
                n_guard = 128;
            otherwise
                error('invalid DRM spectrum occupancy !');
        end

    case 4
        switch(occupancy)
            case 3
                n_guard = 88;
            case 5
                n_guard = 176;
            otherwise
                error('invalid DRM spectrum occupancy !');
        end
    otherwise
        error('invalid DRM robustness mode !');
end