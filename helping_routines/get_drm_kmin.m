%**************************************************************************
% Copyright:     (c) Univerity of Applied Science Rosenheim
% MODULE:        minimum used carrier index
% PROJECT:       DRM Project
% ABBREVIATION:	 drm
% COMPILER:      MATLAB 2006a
% LANGUAGE:      Matlab Interpreter.
% AUTHOR:        Prof. Dr. Markus Stichler
% HISTORY:		 Initial version 11.05.2007
%**************************************************************************
%   returns minimum used carrier index for 
%   DRM Robustness Mode MODE and spectrum occupancy OCCUPANCY
%   References:
%     [1] ETSI: Digital Radio Mondiale System Specification
%**************************************************************************
% function call:
%------------------------
% kmax = get_drm_kmin(mode,occupancy);
%
% Output:
%------------------------
% kmin: minimum used carrier index
%
% Input:
%------------------------
% MODE       : DRM Robustness Mode (A = 1, B = 2, C = 3, D = 4)
% OCCUPANCY  : spectrum occupancy 0 - 5
%**************************************************************************
function kmin = get_drm_kmin(mode,occupancy)

switch(mode)
    case 1
        switch(occupancy)
            case 0
                kmin = 2;
            case 1
                kmin = 2;
            case 2
                kmin = -102;
            case 3
                kmin = -114;
            case 4
                kmin = -98;
            case 5
                kmin = -110;
            otherwise
                error('invalid DRM spectrum occupancy !');
        end
    case 2
        switch(occupancy)
            case 0
                kmin = 1;
            case 1
                kmin = 1;
            case 2
                kmin = -91;
            case 3
                kmin = -103;
            case 4
                kmin = -87;
            case 5
                kmin = -99;
            otherwise
                error('invalid DRM spectrum occupancy !');
        end
    case 3
        switch(occupancy)
            case 3
                kmin = -69;
            case 5
                kmin = -67;
            otherwise
                error('invalid DRM spectrum occupancy !');
        end
    case 4
        switch(occupancy)
            case 3
                kmin = -44;
            case 5
                kmin = -43;
            otherwise
                error('invalid DRM spectrum occupancy !');
        end
    otherwise
        error('invalid DRM robustness mode !');

end
