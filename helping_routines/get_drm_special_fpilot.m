%**************************************************************************
% Copyright:     (c) Univerity of Applied Science Rosenheim
% MODULE:        position of special frequency pilots for DRM Robustness Mode D
% PROJECT:       DRM Project
% ABBREVIATION:	 drm
% COMPILER:      MATLAB 2006a
% LANGUAGE:      Matlab Interpreter.
% AUTHOR:        Prof. Dr. Markus Stichler
% HISTORY:		 Initial version 11.05.2007
%**************************************************************************
%   returns position of special frequency pilots for DRM Robustness Mode D
%   References:
%     [1] ETSI: Digital Radio Mondiale System Specification
%**************************************************************************
% function call:
%------------------------
% special_fpilot = get_drm_special_fpilot;
%
% Output:
%------------------------
% special_fpilot: position of special frequency pilots for DRM Robustness Mode D
%
% Input:
%------------------------
% MODE       : DRM Robustness Mode (A = 1, B = 2, C = 3, D = 4)
% OCCUPANCY  : spectrum occupancy 0 - 5
%**************************************************************************


function special_fpilot = get_drm_special_fpilot

special_fpilot = [7 21];

