function visualize_qo100_link(stSat)
% VISUALIZE_QO100_LINK Visualizes the QO-100 satellite link
%   This function creates visualization for the QO-100 satellite link
%   including orbit, coverage, and signal path
%
%   Parameters:
%   stSat - Structure with satellite parameters

% Make sure we have a scenario
if ~isfield(stSat, 'scenario') || ~isa(stSat.scenario, 'satelliteScenario')
    error('Satellite scenario not found in stSat structure');
end

% Create a viewer for the satellite scenario
v = satelliteScenarioViewer(stSat.scenario, "ShowDetails", true);

% Show orbit trail
showOrbit(v, 'on');

% Add signal path visualization
gs = findobj(stSat.scenario.GroundStations, 'Name', 'Rosenheim');
sat = stSat.scenario.Satellites(1);

% Calculate Earth coverage from QO-100
coverage = conicalSensor(sat, "Name", "NarrowbandTransponder", ...
                         "MountingAngles", [0,0,0], ...
                         "FieldOfView", 17, ...  % Approximately covers Europe, Middle East, Africa
                         "MaxRange", 45000);     % km

% Add coverage visualization
coverageColor = [0.3, 0.6, 0.9, 0.2]; % Light blue with transparency
show(v, coverage, "FaceColor", coverageColor);

% Calculate access times for entire duration
ac = access(gs, sat);
[accessTimes, isAccessible] = accessStatus(ac);

if any(isAccessible)
    fprintf('QO-100 is accessible from ground station\n');
    fprintf('Access intervals:\n');
    for i = 1:size(accessTimes, 1)
        fprintf('  Start: %s, End: %s\n', accessTimes(i,1), accessTimes(i,2));
    end
else
    fprintf('Warning: QO-100 is not accessible from ground station\n');
end

% Draw line of sight from ground station to satellite
showLineOfSight(v, gs, sat, 'Color', 'green', 'LineWidth', 2);

% Show additional information
title(v.Axes, sprintf('QO-100 Satellite Link (%.2f° E)', 25.8));
legend(v.Axes, 'show');

% Create a separate figure for link budget visualization
figure;
tiledlayout(2,2);

% SNR over distance
nexttile;
distances = linspace(35000, 40000, 100);
snrs = zeros(size(distances));
for i = 1:length(distances)
    path_loss = fspl(distances(i)*1000, stSat.uplink_freq);
    snrs(i) = stSat.EIRP - path_loss + stSat.G_T - stSat.k - 10*log10(stSat.bandwidth);
end
plot(distances, snrs, 'LineWidth', 2);
hold on;
plot(stSat.range, stSat.SNR, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
xlabel('Distance to Satellite (km)');
ylabel('SNR (dB)');
title('Link SNR vs. Distance');
grid on;

% Antenna pattern
nexttile;
antenna = dish(stSat.antenna_diameter, stSat.uplink_freq, "Efficiency", stSat.antenna_efficiency);
pattern(antenna, stSat.uplink_freq, -90:90, 0, 'CoordinateSystem', 'polar');
title(sprintf('Dish Antenna Pattern (%.1f m)', stSat.antenna_diameter));

% Link budget breakdown
nexttile;
budget = [stSat.tx_power_dBm, stSat.antenna_gain, -stSat.cable_loss, ...
          -stSat.pointing_loss, -stSat.path_loss, stSat.G_T, -10*log10(stSat.bandwidth)];
labels = {'TX Power', 'Ant Gain', 'Cable Loss', 'Pointing Loss', ...
          'Path Loss', 'G/T', 'BW'};
bar(budget, 'FaceColor', 'flat');
xticklabels(labels);
xtickangle(45);
ylabel('dB');
title('Link Budget Components');
grid on;

% Doppler shift over time
nexttile;
time = stSat.scenario.StartTime:minutes(10):stSat.scenario.StopTime;
doppler = zeros(size(time));

for i = 1:length(time)
    states = states(sat, time(i), "CoordinateFrame", "ECEF");
    pos = states(1:3);
    vel = states(4:6);
    [az, el, r] = lookangles(gs, sat, time(i));
    doppler(i) = satellite.doppler(stSat.uplink_freq, vel, az, el);
end

plot(time, doppler, 'LineWidth', 2);
xlabel('Time');
ylabel('Doppler Shift (Hz)');
title('Doppler Shift Over Time');
grid on;

sgtitle('QO-100 Link Analysis');
end
