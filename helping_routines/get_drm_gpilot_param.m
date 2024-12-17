%**************************************************************************
% Copyright:     (c) Univerity of Applied Science Rosenheim
% MODULE:        gain pilot calculation parameters
% PROJECT:       DRM Project
% ABBREVIATION:	 drm
% COMPILER:      MATLAB 2006a
% LANGUAGE:      Matlab Interpreter.
% AUTHOR:        Prof. Dr. Markus Stichler
% HISTORY:		 Initial version 11.05.2007
%**************************************************************************
%   returns gain pilot calculation parameters for DRM Robustness Mode MODE
%   References:
%     [1] ETSI: Digital Radio Mondiale System Specification
%**************************************************************************
% function call:
%------------------------
% fpilot_position = get_drm_fpilot_position(mode);
%
% Output:
%------------------------
% [X,Y,K0,W1024,Z256,Q1024]: gain pilot calculation parameters
%
% Input:
%------------------------
% MODE       : DRM Robustness Mode (A = 1, B = 2, C = 3, D = 4)
%**************************************************************************

function [x,y,k0,W1024,Z256,Q1024] = get_drm_gpilot_param(mode)


switch(mode)
    case 1
        x = 4;
        y = 5;
        k0 = 2;
        W1024 = [
            228 341 455;
            455 569 683;
            683 796 910;
            910   0 114;
            114 228 341
            ];
        Z256 = [
              0  81 248;
             18 106 106;
            122 116  31;
            129 129  39;
             33  32 111
            ];
        Q1024 = 36;
    case 2
        x = 2;
        y = 3;
        k0 = 1;
        W1024 = [
            512   0 512   0 512;
              0 512   0 512   0;
            512   0 512   0 512;
            ];
        Z256 = [
              0  57 164  64  12;
            168 255 161 106 118;
             25 232 132 233  38;
            ];
        Q1024 = 12;
    case 3
        x = 2;
        y = 2;
        k0 = 1;
        W1024 = [
            465 372 279 186  93   0 931 838 745 652;
            931 838 745 652 559 465 372 279 186  93
            ];
        Z256 = [
              0  76  29  76   9 190 161 248  33 108;
            179 178  83 253 127 105 101 198 250 145
            ];
        Q1024 = 12;
    case 4
        x = 1;
        y = 3;
        k0 = 1;
        W1024 = [
            366 439 512 585 658 731 805 878;
            731 805 878 951   0  73 146 219;
             73 146 219 293 366 439 512 585
            ];
        Z256 = [
              0 240  17  60 220  38 151 101;
            110   7  78  82 175 150 106  25;
            165   7 252 124 253 177 197 142
            ];
        Q1024 = 14;
    otherwise
        error('invalid DRM robustness mode !');
end

