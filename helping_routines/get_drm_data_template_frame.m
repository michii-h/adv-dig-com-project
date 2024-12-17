%**************************************************************************
% Copyright:     (c) Univerity of Applied Science Rosenheim
% MODULE:        drm frame / data symbols = 1 and pilot symbols = 0
% PROJECT:       DRM Project
% ABBREVIATION:	 drm
% COMPILER:      MATLAB 2006a
% LANGUAGE:      Matlab Interpreter.
% AUTHOR:        Prof. Dr. Markus Stichler
% HISTORY:		 Initial version 11.05.2007
%**************************************************************************
%   returns DRM frame with all data symbols = 1 and pilot symbols = 0 for DRM Robustness Mode
%   MODE and spectrum occupancy OCCUPANCY
%   References:
%     [1] ETSI: Digital Radio Mondiale System Specification
%**************************************************************************
% function call:
%------------------------
% data_template_frame = get_drm_data_template_frame(mode, occupancy);
%
% Output:
%------------------------
% data_template_frame: DRM frame with all data symbols = 1 and pilot
% symbols = 0
%
% Input:
%------------------------
% MODE       : DRM Robustness Mode (A = 1, B = 2, C = 3, D = 4)
% OCCUPANCY  : spectrum occupancy 0 - 5
%**************************************************************************
function data_template_frame = get_drm_data_template_frame(mode, occupancy);

% Piloten laden
Plk = get_drm_pilot_frame(mode,occupancy);

symbols_per_frame = get_drm_symbols_per_frame(mode);
n_useful = get_drm_n_useful(mode,occupancy);

dc_position = get_drm_dc_position(mode,occupancy);
kmin = get_drm_kmin(mode,occupancy) + dc_position;
kmax = get_drm_kmax(mode,occupancy) + dc_position;


data_template_frame = ones(symbols_per_frame, n_useful);

data_template_frame(:,1:kmin-1)=0;
data_template_frame(:,kmax+1:end)=0;
data_template_frame(:,dc_position)=0;
data_template_frame(Plk ~=0) = 0;
