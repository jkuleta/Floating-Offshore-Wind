close all; clear all; clc;

% List of result files and corresponding labels
result_files = {'RESULTS_S4.mat', 'RESULTS_S3.mat','RESULTS_S2.mat', 'RESULTS_S1.mat'};
sea_states = {'Sea State 4','Sea State 3', 'Sea State 2', 'Sea State 1'};
motions = {'Surge [m]', 'Heave [m]', 'Pitch [deg]'};
motion_fields = {'Surge','Heave', 'Pitch'}; % Adjust if field names differ

% Preallocate for plotting
data = struct();

% Load data from each file
for i = 1:length(result_files)
    S = load(result_files{i});
    for j = 1:length(motion_fields)
        data(i).(motion_fields{j}) = S.RESULTS.(motion_fields{j});
    end
    data(i).time = S.RESULTS.Time;
end

% Plot settings
figure('Position', [100, 100, 600, 500]); % Smaller figure
set(gcf, 'Color', 'w');
fontsize = 10;

t = tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

for j = 1:length(motion_fields)
    ax = nexttile;
    hold on
    for i = 1:length(result_files)
        plot(data(i).time, data(i).(motion_fields{j}), 'LineWidth', 1.5)
    end
    ylabel(motions{j}, 'FontSize', fontsize)
    title(motion_fields{j}, 'FontSize', fontsize)
    set(gca, 'FontSize', fontsize)
    grid on
    if j == length(motion_fields)
        xlabel('Time [s]', 'FontSize', fontsize)
    end
end

% Add legend at the bottom
lgd = legend(sea_states, 'Orientation', 'horizontal', 'FontSize', fontsize);
lgd.Layout.Tile = 'south';

% Save the first figure
saveas(gcf, 'TASK 11 - Mild_motions_plot.png');

% --- Plot Power and Thrust ---
power_thrust_fields = {'Power', 'Thrust'};
power_thrust_labels = {'Power [MW]', 'Thrust [MN]'};

% Load data from each file (if not already loaded)
for i = 1:length(result_files)
    S = load(result_files{i});
    for j = 1:length(power_thrust_fields)
        data(i).(power_thrust_fields{j}) = S.RESULTS.(power_thrust_fields{j});
    end
    data(i).time = S.RESULTS.Time;
end

% Create figure with tiledlayout for equal-sized plots
figure('Position', [100, 100, 600, 350]);
set(gcf, 'Color', 'w');
fontsize = 10;

t2 = tiledlayout(2,1, 'TileSpacing', 'compact', 'Padding', 'compact');

for j = 1:length(power_thrust_fields)
    ax = nexttile;
    hold on
    for i = 1:length(result_files)
        % Convert units if needed (assume Power in W, Thrust in N)
        if strcmp(power_thrust_fields{j}, 'Power')
            y = data(i).Power / 1e6; % MW
        else
            y = data(i).Thrust / 1e6; % MN
        end
        plot(data(i).time, y, 'LineWidth', 1.5)
    end
    ylabel(power_thrust_labels{j}, 'FontSize', fontsize)
    title(power_thrust_fields{j}, 'FontSize', fontsize)
    set(gca, 'FontSize', fontsize)
    grid on
    xlim([0 599])
    if j == length(power_thrust_fields)
        xlabel('Time [s]', 'FontSize', fontsize)
    end
end

% Add legend across the bottom
lgd2 = legend(sea_states, 'Orientation', 'horizontal', 'FontSize', fontsize);
lgd2.Layout.Tile = 'south';

% Save the Power/Thrust figure
print(gcf, 'TASK 11 - Mild_power_thrust_plot', '-dpng');
