%**************************************************************************
% Copyright:     (c) Univerity of Applied Science Rosenheim
% MODULE:        maximum used carrier index
% PROJECT:       DRM Project
% ABBREVIATION:	 drm
% COMPILER:      MATLAB 2006a
% LANGUAGE:      Matlab Interpreter.
% AUTHOR:        Prof. Dr. Markus Stichler
% HISTORY:		 Initial version 11.05.2007
%**************************************************************************
%   returns maximum used carrier index for 
%   DRM Robustness Mode MODE and spectrum occupancy OCCUPANCY
%   References:
%     [1] ETSI: Digital Radio Mondiale System Specification
%**************************************************************************
% function call:
%------------------------
% kmax = get_drm_kmax(mode,occupancy);
%
% Output:
%------------------------
% kmax: maximum used carrier index
%
% Input:
%------------------------
% MODE       : DRM Robustness Mode (A = 1, B = 2, C = 3, D = 4)
% OCCUPANCY  : spectrum occupancy 0 - 5
%**************************************************************************

function kmax = get_drm_kmax(mode,occupancy)

switch(mode)
    case 1
        switch(occupancy)
            case 0
                kmax = 102;
            case 1
                kmax = 114;
            case 2
                kmax = 102;
            case 3
                kmax = 114;
            case 4
                kmax = 314;
            case 5
                kmax = 350;
            otherwise
                error('invalid DRM spectrum occupancy !');
        end
    case 2
        switch(occupancy)
            case 0
                kmax = 91;
            case 1
                kmax = 103;
            case 2
                kmax = 91;
            case 3
                kmax = 103;
            case 4
                kmax = 279;
            case 5
                kmax = 311;
            otherwise
                error('invalid DRM spectrum occupancy !');
        end
    case 3
        switch(occupancy)
            case 3
                kmax = 69;
            case 5
                kmax = 213;
            otherwise
                error('invalid DRM spectrum occupancy !');
        end
    case 4
        switch(occupancy)
            case 3
                kmax = 44;
            case 5
                kmax = 135;
            otherwise
                error('invalid DRM spectrum occupancy !');
        end
    otherwise
        error('invalid DRM robustness mode !');
end
