% STARTUP Set Matlab search path for OFDM Workshop.

%==================================================================================
% Init
%==================================================================================

close all;
clear all;
clear functions;
clear classes;
clc;

fprintf('=======================================================================\n');
fprintf('                            OFDM-Workshop                              \n');
fprintf('=======================================================================\n');


%==================================================================================
% Adding paths to MATLAB directory search path
%==================================================================================

% Get current path
ActualPath = pwd;

% Add ORION testbench
fprintf('* Adding DRM to MATLAB directory search path ...\n')
addpath(genpath(ActualPath));

