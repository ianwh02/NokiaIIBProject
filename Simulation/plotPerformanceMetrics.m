% plotPerformanceMetrics.m
% Updated for structured .mat files with simResults.{simThroughput, maxThroughput, bler}

clear; clc; close all;

% === MANUALLY SPECIFY YOUR FILES HERE ===
resultsFolder = './results/';
fileNames = {
    'simResult_HARQ_perfectCSI_perfectestimator_CDL-C_8Tx_4Rx_4Layer_10Hz.mat'
    % Add more filenames here
};

% Define the SNR vector (based on your simulation sweep)
snrVec = -5:1:25;

% Plot formatting
colors = lines(numel(fileNames));
markers = {'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h', 'x'};

% Create figures
fig1 = figure; hold on; grid on;
xlabel('SNR (dB)'); ylabel('Throughput (%)');
title('Throughput vs. SNR');
set(gca, 'FontSize', 12);

fig2 = figure; hold on; grid on;
xlabel('SNR (dB)'); ylabel('BLER');
title('BLER vs. SNR');
set(gca, 'YScale', 'log', 'FontSize', 12);
yline(0.1, '--r', 'Target BLER = 10%', 'LabelHorizontalAlignment', 'left');

% Loop through all specified files
for i = 1:length(fileNames)
    filePath = fullfile(resultsFolder, fileNames{i});
    data = load(filePath);

    % Extract and transpose as needed
    simThroughput = double(data.simResults.simThroughput(:));  % 1xN
    maxThroughput = double(data.simResults.maxThroughput(:));
    bler = double(data.simResults.bler(:));

    % Compute throughput %
    throughputPercent = (simThroughput ./ maxThroughput) * 100;

    % Legend name from file name (remove .mat)
    legendName = erase(fileNames{i}, '.mat');

    % Plot Throughput
    figure(fig1);
    plot(snrVec, throughputPercent, ...
        'LineWidth', 1.8, ...
        'Color', colors(i,:), ...
        'Marker', markers{mod(i-1,numel(markers))+1}, ...
        'DisplayName', legendName);

    % Plot BLER
    figure(fig2);
    semilogy(snrVec, bler, ...
        'LineWidth', 1.8, ...
        'Color', colors(i,:), ...
        'Marker', markers{mod(i-1,numel(markers))+1}, ...
        'DisplayName', legendName);
end

% Final formatting
figure(fig1);
legend('Location', 'best');
ylim([0 100]);

figure(fig2);
legend('Location', 'best');
ylim([1e-3 1]);