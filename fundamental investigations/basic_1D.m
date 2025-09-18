% --- Initialize Structure ---
close all;

savePath = 'G:\Meine Ablage\Studium\Master\Project Modulatory Feedback\FG_RNN\output\feedback\20250820_Negative_B_Cells\high resolution\long C';
if ~exist(savePath, 'dir')
    mkdir(savePath);
end

% --- Parameters ---
iterations = 100;
alphaB = 0.36;
betaB = 1;
gammaB = 7;
zetaB = 0.7;
lambdaB = 5;
xiB = 0;

alphaG = 0.25;  
betaG = 1;
gammaG = 0;
zetaG = 0.9;
muG = 1;

% --- Initialize Structure ---
[B_in_cell, B_out_cell, G_vec, W_in_cell, W_out_cell] = build1DStructure(getStringMap(7), 6, true, fullfile(savePath, "spatial structure.png"));

num_shapes = length(B_in_cell);
nG = length(G_vec);

% Get total number of B cells across all shapes
total_nB = 0;
nB_per_shape = zeros(num_shapes, 1);
for shape_idx = 1:num_shapes
    nB_per_shape(shape_idx) = length(B_in_cell{shape_idx});
    total_nB = total_nB + nB_per_shape(shape_idx);
end

% --- Initialize Histories ---
% Store histories for each shape separately
B_in_history_cell = cell(num_shapes, 1);
B_out_history_cell = cell(num_shapes, 1);
B_in_prev_cell = cell(num_shapes, 1);
B_out_prev_cell = cell(num_shapes, 1);

for shape_idx = 1:num_shapes
    B_in_history_cell{shape_idx} = zeros(nB_per_shape(shape_idx), iterations+1);
    B_out_history_cell{shape_idx} = zeros(nB_per_shape(shape_idx), iterations+1);
    B_in_prev_cell{shape_idx} = zeros(nB_per_shape(shape_idx), 1);
    B_out_prev_cell{shape_idx} = zeros(nB_per_shape(shape_idx), 1);
    
    % Initial states
    B_in_history_cell{shape_idx}(:,1) = B_in_prev_cell{shape_idx};
    B_out_history_cell{shape_idx}(:,1) = B_out_prev_cell{shape_idx};
end

G_history = zeros(nG, iterations+1);
G_history(:,1) = G_vec;

% --- Main Simulation Loop ---
for t = 1:iterations
    % Initialize total G input
    total_G_input = zeros(nG, 1);
    total_G_inh = zeros(nG, 1);
    
    % Process each shape
    for shape_idx = 1:num_shapes
        % Feedback from G to B for this shape
        FB_in = W_in_cell{shape_idx}' * G_vec;
        FB_out = W_out_cell{shape_idx}' * G_vec;

        % Update B_in and B_out for this shape
        B_in_current = betaB * (1 + lambdaB * FB_in) ./ ( ...
            alphaB + gammaB * B_out_prev_cell{shape_idx} + zetaB * (1 + lambdaB * FB_in));
        
        B_out_current = betaB * (1 + lambdaB * FB_out) ./ ( ...
            alphaB + gammaB * B_in_prev_cell{shape_idx} + zetaB * (1 + lambdaB * FB_out));

        % Store updated B activities
        B_in_cell{shape_idx} = B_in_current;
        B_out_cell{shape_idx} = B_out_current;
        
        % Accumulate input from this shape to G
        shape_G_input = W_in_cell{shape_idx} * B_in_current + W_out_cell{shape_idx} * B_out_current;
        total_G_input = total_G_input + shape_G_input;

        shape_G_inh = W_out_cell{shape_idx} * B_in_current + W_in_cell{shape_idx} * B_out_current;
        total_G_inh = total_G_inh + shape_G_inh;
        
        % Store history for this shape
        B_in_history_cell{shape_idx}(:,t+1) = B_in_current;
        B_out_history_cell{shape_idx}(:,t+1) = B_out_current;
        
        % Update previous values for this shape
        B_in_prev_cell{shape_idx} = B_in_current;
        B_out_prev_cell{shape_idx} = B_out_current;
    end
    
    % Update G cells with total input from all shapes
    G_norm = alphaG + gammaG * total_G_input + zetaG * total_G_input + muG * total_G_inh;
    G_vec = betaG * total_G_input ./ G_norm;
    
    % Store G history
    G_history(:,t+1) = G_vec;
end

%% --- Plot Final Activities ---
figure('Name', 'Activities at Different Iterations');

% Define which iterations to plot (every 5th iteration plus the final one)
n = 7;
iterations_to_show = [1:round(iterations/n):iterations, iterations+1];
if length(iterations_to_show) > 1 && iterations_to_show(end-1) == iterations_to_show(end)
    iterations_to_show(end-1) = []; % Remove duplicate if iterations is multiple of 5
end

% Create colormap for different iterations
colors = copper(length(iterations_to_show));

% Create subplots: one row for each shape, plus one for G cells
total_subplots = num_shapes + 1;

% Store handles for legend
legend_handles = [];
legend_labels = {};

gap = 0.1;  % Space between subplots
available_height = 0.8 - (total_subplots - 1) * gap;
subplot_height = available_height / total_subplots;

% Plot B cells for each shape separately
for shape_idx = 1:num_shapes
    left = 0.1;
    width = 0.7;
    height = subplot_height;
    bottom = 0.1 + (total_subplots - shape_idx) * (subplot_height + gap);
    
    subplot('Position', [left, bottom, width, height]);
    hold on;
    for i = 1:length(iterations_to_show)
        iter = iterations_to_show(i);
        h1 = plot(B_in_history_cell{shape_idx}(:, iter), '.-', 'Color', colors(i,:), ...
            'LineWidth', 1.5);
        h2 = plot(B_out_history_cell{shape_idx}(:, iter), '.--', 'Color', colors(i,:), ...
            'LineWidth', 1.5);
        
        % Store handles for legend from first shape only
        if shape_idx == 1
            legend_handles = [legend_handles, h1, h2];
            legend_labels = [legend_labels, sprintf('B_{in} (iter %d)', iter-1), sprintf('B_{out} (iter %d)', iter-1)];
        end
    end
    title(sprintf('Shape %d - B cell activities at different iterations (nB = %d)', shape_idx, nB_per_shape(shape_idx)));
    xlabel('B cell index within shape');
    ylabel('Activity');
    grid on;
end

% Plot G cells
left = 0.1;
width = 0.7;
height = subplot_height;
bottom = 0.1;  % Bottom subplot

subplot('Position', [left, bottom, width, height]);
hold on;
for i = 1:length(iterations_to_show)
    iter = iterations_to_show(i);
    h3 = plot(G_history(:, iter), '.-', 'Color', colors(i,:), ...
        'LineWidth', 1.5);
    
    % Add G handles to legend
    legend_handles = [legend_handles, h3];
    legend_labels = [legend_labels, sprintf('G (iter %d)', iter-1)];
end
title('G cell activities at different iterations');
xlabel('G cell index');
ylabel('Activity');
grid on;

% Add local legends to B cell subplots
for shape_idx = 1:num_shapes
    bottom_pos = 0.1 + (total_subplots - shape_idx) * (subplot_height + gap);
    subplot('Position', [0.1, bottom_pos, 0.7, subplot_height]);
    h_in = plot(NaN, NaN, '-', 'Color', 'k', 'LineWidth', 1.5);
    h_out = plot(NaN, NaN, '--', 'Color', 'k', 'LineWidth', 1.5);
    legend([h_in, h_out], {'B_{in}', 'B_{out}'}, 'Location', 'best', 'FontSize', 8);
end

% Create global colorbar for iterations
colormap(copper(length(iterations_to_show)));
c = colorbar('eastoutside');
c.Label.String = 'Iteration';
c.Ticks = linspace(0, 1, length(iterations_to_show));
c.TickLabels = arrayfun(@(x) sprintf('%d', x-1), iterations_to_show, 'UniformOutput', false);
c.Position = [0.85, 0.1, 0.03, 0.8];

set(gcf, 'Position', [100, 100, 1200, 300*total_subplots]);
print(fullfile(savePath, 'total_activities.png'), '-dpng', '-r300');

%% --- Plot Activities by Shape ---
figure('Name', 'Activities by Shape');
iterations_axis = 1:iterations+1;

% Plot each shape separately
for shape_idx = 1:num_shapes
    subplot(num_shapes, 1, shape_idx);
    hold on;
    
    % Plot final iteration for this shape
    plot(B_in_history_cell{shape_idx}(:, end), 'r.-', 'LineWidth', 1.5, 'MarkerSize', 8, 'DisplayName', 'B_{in}');
    plot(B_out_history_cell{shape_idx}(:, end), 'b.-', 'LineWidth', 1.5, 'MarkerSize', 8, 'DisplayName', 'B_{out}');
    
    title(sprintf('Shape %d - Final Activities (nB = %d)', shape_idx, nB_per_shape(shape_idx)));
    xlabel('B cell index within shape');
    ylabel('Activity');
    legend('Location', 'northeast');
    grid on;
end

set(gcf, 'Position', [100, 100, 1000, 200*num_shapes]);
print(fullfile(savePath, 'activities_by_shape.png'), '-dpng', '-r300');

%% --- Individual B Cell Traces ---
% Automatically select some B cells from each shape for detailed plotting
B_indices_to_plot_per_shape = cell(num_shapes, 1);
for shape_idx = 1:num_shapes
    nB_shape = nB_per_shape(shape_idx);
    if nB_shape <= 8
        B_indices_to_plot_per_shape{shape_idx} = 1:nB_shape;
    else
        % Select evenly spaced indices
        indices = round(linspace(1, nB_shape, 8));
        B_indices_to_plot_per_shape{shape_idx} = unique(indices);
    end
end
B_indices_to_plot_per_shape{1} = [100, 130, 140, 150, 165, 168, 171, 175];

% Plot B cell traces for each shape
maxPerRow = 4;
for shape_idx = 1:num_shapes
    indices_to_plot = B_indices_to_plot_per_shape{shape_idx};
    numB = length(indices_to_plot);
    numBRows = ceil(numB / maxPerRow);
    
    figure('Name', sprintf('Shape %d B Cell Activities', shape_idx));
    
    % Find global limits for this shape
    B_max = max([B_in_history_cell{shape_idx}(:); B_out_history_cell{shape_idx}(:)]);
    B_ylim = [0, B_max + 0.05*B_max];
    
    for idx = 1:numB
        b_idx = indices_to_plot(idx);
        B_in_trace = B_in_history_cell{shape_idx}(b_idx, :);
        B_out_trace = B_out_history_cell{shape_idx}(b_idx, :);
        
        subplot(numBRows, min(maxPerRow, numB), idx);
        plot(iterations_axis, B_in_trace, 'r-', 'LineWidth', 1.5); hold on;
        plot(iterations_axis, B_out_trace, 'b-', 'LineWidth', 1.5);
        title(sprintf('Shape %d, B Cell %d', shape_idx, b_idx));
        xlabel('Iter'); ylabel('Activity');
        legend('B_{in}','B_{out}','Location','best');
        grid on;
        ylim(B_ylim);
    end
    sgtitle(sprintf('Shape %d B Cell Activities Over Iterations', shape_idx));
    set(gcf, "Position", [100, 100, 1200, 300*numBRows]);
    print(fullfile(savePath, sprintf('B_cell_activities_shape_%d.png', shape_idx)), '-dpng', '-r300');
end

%% --- G Cells Activity ---
% Select some G cells for detailed plotting
if nG <= 8
    G_indices_to_plot = 1:nG;
else
    G_indices_to_plot = round(linspace(1, nG, 8));
    G_indices_to_plot = unique(G_indices_to_plot);
end
G_indices_to_plot = [607, 617, 642, 944, 949, 960, 970, 979];

numG = length(G_indices_to_plot);
numGRows = ceil(numG / maxPerRow);

% Find global limits for G cells
G_max = max(G_history(:));
G_ylim = [0, G_max + 0.05*G_max];

figure('Name', 'G Cell Activities');
for idx = 1:numG
    g_idx = G_indices_to_plot(idx);
    G_trace = G_history(g_idx, :);
    subplot(numGRows, min(maxPerRow, numG), idx);
    plot(iterations_axis, G_trace, 'k-', 'LineWidth', 1.5);
    title(sprintf('G Cell Index %d', g_idx));
    xlabel('Iter'); ylabel('Activity');
    legend('G','Location','best');
    grid on;
    ylim(G_ylim);
end
sgtitle('G Cell Activities Over Iterations');
set(gcf, "Position", [100, 100, 1200, 300*numGRows]);
print(fullfile(savePath, 'G_cell_activities.png'), '-dpng', '-r300');

%% Parameters
% --- Parameter Figure ---
figure('Name', 'Parameters', 'Position', [100, 100, 450, 700]);

% Create parameter string with proper formatting
paramText = sprintf([...
    '\\fontsize{16}\\bfModel Parameters\n\n' ...
    '\\fontsize{12}\\bfStructure:\\rm\n' ...
    '    Number of shapes = %d\n' ...
    '    Total B cells = %d\n' ...
    '    Total G cells = %d\n\n' ...
    '\\fontsize{12}\\bfB Cell Parameters:\\rm\n' ...
    '    \\alpha_B = %.2f\n' ...
    '    \\beta_B = %.2f\n' ...
    '    \\gamma_B = %.2f\n' ...
    '    \\zeta_B = %.2f\n' ...
    '    \\lambda_B = %.2f\n' ...
    '    \\xi_B = %.2f\n\n' ...
    '\\fontsize{12}\\bfG Cell Parameters:\\rm\n' ...
    '    \\alpha_G = %.2f\n' ...
    '    \\beta_G = %.2f\n' ...
    '    \\gamma_G = %.2f\n' ...
    '    \\zeta_G = %.2f\n' ...
    '    \\mu_G = %.2f\n\n' ...
    '\\fontsize{12}\\bfSimulation:\\rm\n' ...
    '    Iterations = %d'], ...
    num_shapes, total_nB, nG, ...
    alphaB, betaB, gammaB, zetaB, lambdaB, xiB, ...
    alphaG, betaG, gammaG, zetaG, muG, ...
    iterations);

% Create axes that fills the figure
ax = axes('Position', [0.05 0.05 0.9 0.9], 'Visible', 'off');
xlim([0 1]);
ylim([0 1]);

% Add background rectangle with padding
rectangle('Position', [0.05 0.05 0.9 0.9], ...
    'FaceColor', [0.95 0.95 0.95], ...
    'EdgeColor', 'black', ...
    'LineWidth', 2);

% Add text with proper positioning
text(0.5, 0.95, paramText, ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'top', ...
    'FontName', 'FixedWidth', ...
    'FontSize', 11, ...
    'Interpreter', 'tex');

axis off;
print(fullfile(savePath, 'parameters.png'), '-dpng', '-r300');