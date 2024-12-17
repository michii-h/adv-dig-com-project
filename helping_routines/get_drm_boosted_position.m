%**************************************************************************
% Copyright:     (c) Univerity of Applied Science Rosenheim
% MODULE:        boosted gain pilot positions
% PROJECT:       DRM Project
% ABBREVIATION:	 drm
% COMPILER:      MATLAB 2006a
% LANGUAGE:      Matlab Interpreter.
% AUTHOR:        Prof. Dr. Markus Stichler
% HISTORY:		 Initial version 11.05.2007
%**************************************************************************
% returns boosted gain pilot positions for DRM Robustness Mode MODE and
% spectrum occupancy OCCUPANCY
%   References:
%     [1] ETSI: Digital Radio Mondiale System Specification
%**************************************************************************
% function call:
%------------------------
% boosted_position =  get_drm_boosted_position(mode,occupancy);
%
% Output:
%------------------------
% boosted_position: boosted gain pilot positions
%
% Input:
%------------------------
% MODE       : DRM Robustness Mode (A = 1, B = 2, C = 3, D = 4)
% OCCUPANCY  : spectrum occupancy 0 - 5
%**************************************************************************
function boosted_position =  get_drm_boosted_position(mode,occupancy)

switch(mode)
    case 1
        switch(occupancy)
            case 0
                boosted_position = [2 6 98 102];
            case 1
                boosted_position = [2 6 110 114];
            case 2
                boosted_position = [-102 -98 98 102];
            case 3
                boosted_position = [-114 -110 110 114];
            case 4
                boosted_position = [-98 -94 310 314];
            case 5
                boosted_position = [-110 -106 346 350];
            otherwise
                error('invalid DRM spectrum occupancy !');
        end
    case 2
        switch(occupancy)
            case 0
                boosted_position = [1 3 89 91];
            case 1
                boosted_position = [1 3 101 103];
            case 2
                boosted_position = [-91 -89 89 91];
            case 3
                boosted_position = [-103 -101 101 103];
            case 4
                boosted_position = [-87 -85 277 279];
            case 5
                boosted_position = [-99 -97 309 311];
            otherwise
                error('invalid DRM spectrum occupancy !');
        end
    case 3
        switch(occupancy)
            case 3
                boosted_position = [-69 -67 67 69];
            case 5
                boosted_position = [-67 -65 211 213];
            otherwise
                error('invalid DRM spectrum occupancy !');
        end
    case 4
        switch(occupancy)
            case 3
                boosted_position = [-44 -43 43 44];
            case 5
                boosted_position = [-43 -42 134 135];
            otherwise
                error('invalid DRM spectrum occupancy !');
        end
    otherwise
        error('invalid DRM robustness mode !');
end
