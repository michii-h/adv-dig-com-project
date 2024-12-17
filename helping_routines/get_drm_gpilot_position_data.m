%**************************************************************************
% Copyright:     (c) Univerity of Applied Science Rosenheim
% MODULE:        gain pilot positions and data
% PROJECT:       DRM Project
% ABBREVIATION:	 drm
% COMPILER:      MATLAB 2006a
% LANGUAGE:      Matlab Interpreter.
% AUTHOR:        Prof. Dr. Markus Stichler
% HISTORY:		 Initial version 11.05.2007
%**************************************************************************
%   returns gain pilot positions and data for DRM Robustness Mode MODE and
%   spectrum occupancy OCCUPANCY
%   References:
%     [1] ETSI: Digital Radio Mondiale System Specification
%**************************************************************************
% function call:
%------------------------
% [gpilot_position,gpilot_data] = get_drm_gpilot_position_data(mode,occupancy);
%
% Output:
%------------------------
% gpilot_position: gain pilot positions
% gpilot_data: gain pilot data
%
% Input:
%------------------------
% MODE       : DRM Robustness Mode (A = 1, B = 2, C = 3, D = 4)
% OCCUPANCY  : spectrum occupancy 0 - 5
%**************************************************************************

function [gpilot_position,gpilot_data] = get_drm_gpilot_position_data(mode,occupancy)

[x,y,k0,W1024,Z256,Q1024] = get_drm_gpilot_param(mode);
symbols_per_frame = get_drm_symbols_per_frame(mode);
kmin = get_drm_kmin(mode,occupancy);
kmax = get_drm_kmax(mode,occupancy);

magnitude = sqrt(2);
increment = pi / 512;

for s = 1:symbols_per_frame
    n = mod((s-1),y);
    m = floor((s-1)/y);
    pilot_count = 0;
    for k = kmin:kmax
        p = (k - k0 - n*x)/(x*y);
        if p == floor(p) % Falls Gain-Pilot Frequenz, entsprechende Phase berechnen
            pilot_count = pilot_count + 1;
            gpilot_position(s,pilot_count) = k;
            gpilot_phase = mod(4*Z256(n+1,m+1) + p*W1024(n+1,m+1) + p^2*s*Q1024,1024);
            gpilot_data(s,pilot_count) = magnitude * exp(j * gpilot_phase * increment);
        end
    end
end

