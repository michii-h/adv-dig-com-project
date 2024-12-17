%**************************************************************************
% Copyright:     (c) Univerity of Applied Science Rosenheim
% MODULE:        timing pilot data
% PROJECT:       DRM Project
% ABBREVIATION:	 drm
% COMPILER:      MATLAB 2006a
% LANGUAGE:      Matlab Interpreter.
% AUTHOR:        Prof. Dr. Markus Stichler
% HISTORY:		 Initial version 11.05.2007
%**************************************************************************
%  returns timing pilot data for DRM Robustness Mode MODE
%   References:
%     [1] ETSI: Digital Radio Mondiale System Specification
%**************************************************************************
% function call:
%------------------------
% tpilot_data = get_drm_tpilot_data(mode);
%
% Output:
%------------------------
% tpilot_data: timing pilot data
%
% Input:
%------------------------
% MODE       : DRM Robustness Mode (A = 1, B = 2, C = 3, D = 4)
%**************************************************************************


function tpilot_data = get_drm_tpilot_data(mode)


magnitude = sqrt(2);
increment = pi / 512;

switch(mode)
    case 1
        tpilot_data = [
            magnitude * exp(j * 973 * increment)
            magnitude * exp(j * 717 * increment)
            magnitude * exp(j * 264 * increment)
            magnitude * exp(j * 357 * increment)
            magnitude * exp(j * 357 * increment)
            magnitude * exp(j * 952 * increment)
            magnitude * exp(j * 440 * increment)
            magnitude * exp(j * 856 * increment)
            magnitude * exp(j *  88 * increment)
            magnitude * exp(j *  88 * increment)
            magnitude * exp(j *  68 * increment)
            magnitude * exp(j * 836 * increment)
            magnitude * exp(j * 836 * increment)
            magnitude * exp(j *1008 * increment)
            magnitude * exp(j *1008 * increment)
            magnitude * exp(j * 752 * increment)
            magnitude * exp(j * 215 * increment)
            magnitude * exp(j * 727 * increment)
            ];
    case 2
        tpilot_data = [
            magnitude * exp(j * 304 * increment)
            magnitude * exp(j * 108 * increment)
            magnitude * exp(j * 620 * increment)
            magnitude * exp(j * 192 * increment)
            magnitude * exp(j * 704 * increment)
            magnitude * exp(j *  44 * increment)
            magnitude * exp(j * 432 * increment)
            magnitude * exp(j * 588 * increment)
            magnitude * exp(j * 844 * increment)
            magnitude * exp(j * 651 * increment)
            magnitude * exp(j * 651 * increment)
            magnitude * exp(j * 460 * increment)
            magnitude * exp(j * 460 * increment)
            magnitude * exp(j * 944 * increment)
            magnitude * exp(j * 940 * increment)
            magnitude * exp(j * 428 * increment)
            ];
    case 3
        tpilot_data = [
            magnitude * exp(j * 722 * increment)
            magnitude * exp(j * 466 * increment)
            magnitude * exp(j * 214 * increment)
            magnitude * exp(j * 479 * increment)
            magnitude * exp(j * 516 * increment)
            magnitude * exp(j * 260 * increment)
            magnitude * exp(j * 577 * increment)
            magnitude * exp(j * 662 * increment)
            magnitude * exp(j *   3 * increment)
            magnitude * exp(j * 771 * increment)
            magnitude * exp(j * 392 * increment)
            magnitude * exp(j *  37 * increment)
            magnitude * exp(j *  37 * increment)
            magnitude * exp(j * 474 * increment)
            magnitude * exp(j * 242 * increment)
            magnitude * exp(j * 754 * increment)
            ];
    case 4
        tpilot_data = [
            magnitude * exp(j * 636 * increment)
            magnitude * exp(j * 124 * increment)
            magnitude * exp(j * 788 * increment)
            magnitude * exp(j * 200 * increment)
            magnitude * exp(j * 688 * increment)
            magnitude * exp(j * 152 * increment)
            magnitude * exp(j * 920 * increment)
            magnitude * exp(j * 920 * increment)
            magnitude * exp(j * 644 * increment)
            magnitude * exp(j * 388 * increment)
            magnitude * exp(j * 652 * increment)
            magnitude * exp(j * 176 * increment)
            magnitude * exp(j * 176 * increment)
            magnitude * exp(j * 752 * increment)
            magnitude * exp(j * 496 * increment)
            magnitude * exp(j * 432 * increment)
            magnitude * exp(j * 964 * increment)
            magnitude * exp(j * 452 * increment)
            ];
    otherwise
        error('invalid DRM robustness mode !');
end
