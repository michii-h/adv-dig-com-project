%**************************************************************************
% Copyright:     (c) Univerity of Applied Science Rosenheim
% MODULE:        length of useful symbol
% PROJECT:       DRM Project
% ABBREVIATION:	 drm
% COMPILER:      MATLAB 2006a
% LANGUAGE:      Matlab Interpreter.
% AUTHOR:        Prof. Dr. Markus Stichler
% HISTORY:		 Initial version 11.05.2007
%**************************************************************************
%   returns length of useful symbol
%   part for DRM Robustness Mode MODE and spectrum occupancy OCCUPANCY
%   References:
%     [1] ETSI: Digital Radio Mondiale System Specification
%**************************************************************************
% function call:
%------------------------
% n_useful = get_drm_n_useful(mode,occupancy);
%
% Output:
%------------------------
% n_useful: length of useful symbol
%
% Input:
%------------------------
% MODE       : DRM Robustness Mode (A = 1, B = 2, C = 3, D = 4)
% OCCUPANCY  : spectrum occupancy 0 - 5
%**************************************************************************

function n_useful = get_drm_n_useful(mode,occupancy)

switch(mode)
    case 1
        switch(occupancy)
            case {0,1,2,3}
                n_useful = 288;
            case {4,5}
                n_useful = 576;
            otherwise
                error('invalid DRM spectrum occupancy !');
        end
    case 2
        switch(occupancy)
            case {0,1,2,3}
                n_useful = 256;
            case {4,5}
                n_useful = 512;
            otherwise
                error('invalid DRM spectrum occupancy !');
        end
    case 3
        switch(occupancy)
            case 3
                n_useful = 176;
            case 5
                n_useful = 352;
            otherwise
                error('invalid DRM spectrum occupancy !');
        end
    case 4
        switch(occupancy)
            case 3
                n_useful = 112;
            case 5
                n_useful = 224;
            otherwise
                error('invalid DRM spectrum occupancy !');
        end
    otherwise
        error('invalid DRM robustness mode !');
end
