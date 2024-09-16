load("test_feedback.mat");

% Define parameters
params.bPrs.exp_decay = 3;        % alpha
params.bPrs.FF.inhibition = 1;    % delta
params.bPrs.FF.spatial_neighborhood = gaussianFilter2D(13, 13, 1, 1);

params.bPrs.FB.scale = 1;           % lambda
params.bPrs.FB.offset = 0;          % T_offset
params.bPrs.FB.inhibition = 1;    % gamma
params.bPrs.FB.spatial_neighborhood = gaussianFilter2D(13, 13, 1, 1);

[bPyr1_1_FB, bPyr1_2_FB, bPyr2_1_FB, bPyr2_2_FB, cPyr_FB] = grouping_feedback(gPyr1, gPyr2, bPyr1_1, bPyr1_2, bPyr2_1, bPyr2_2, cPyr, params);

for ori = 1:8
    all_data = [bPyr1_1(1).orientation(ori).data(:); bPyr1_2(1).orientation(ori).data(:); ...
                bPyr2_1(1).orientation(ori).data(:); bPyr2_2(1).orientation(ori).data(:); ...
                bPyr1_1_FB(1).orientation(ori).data(:); bPyr1_2_FB(1).orientation(ori).data(:); ...
                bPyr2_1_FB(1).orientation(ori).data(:); bPyr2_2_FB(1).orientation(ori).data(:)];
    global_min = min(all_data);
    global_max = max(all_data);

    figure;
    sgtitle(['Orientation ', num2str(ori)]); % Title for the figure

    % Create a 4x4 tiled layout
    t = tiledlayout(4, 4, 'TileSpacing', 'tight', 'Padding', 'tight'); % Adjust spacing

    % Plot the G-cell activations for Light (gPyr1)
    nexttile([2, 2]); % Make G-cell activity span over two columns
    imagesc(gPyr1(1).data); % G-cell activations L 1_1
    title('G-Cells L');
    axis off; colormap(gray); colorbar;

    % First row for Light (L) B-cells
    nexttile;
    imagesc(bPyr1_1(1).orientation(ori).data, [global_min global_max]); % B-cells L theta
    title('B-Cells L theta');
    axis off; colormap(gray); colorbar;

    nexttile;
    imagesc(bPyr1_2(1).orientation(ori).data, [global_min global_max]); % B-cells L theta+pi
    title('B-Cells L theta+pi');
    axis off; colormap(gray); colorbar;

    % Second row for Light (L) FB-cells
    nexttile;
    imagesc(bPyr1_1_FB(1).orientation(ori).data, [global_min global_max]); % FB-cells L theta
    title('FB B-Cells L theta');
    axis off; colormap(gray); colorbar;

    nexttile;
    imagesc(bPyr1_2_FB(1).orientation(ori).data, [global_min global_max]); % FB-cells L theta+pi
    title('FB B-Cells L theta+pi');
    axis off; colormap(gray); colorbar;

    % Plot the G-cell activations for Dark (gPyr2)
    nexttile([2, 2]);
    imagesc(gPyr2(1).data); % G-cell activations D 2_1
    title('G-Cells D');
    axis off; colormap(gray); colorbar;

    % Third row for Dark (D) B-cells
    nexttile;
    imagesc(bPyr2_1(1).orientation(ori).data, [global_min global_max]); % B-cells D theta
    title('B-Cells D theta');
    axis off; colormap(gray); colorbar;

    nexttile;
    imagesc(bPyr2_2(1).orientation(ori).data, [global_min global_max]); % B-cells D theta+pi
    title('B-Cells D theta+pi');
    axis off; colormap(gray); colorbar;

    % Fourth row for Dark (D) FB-cells
    nexttile;
    imagesc(bPyr2_1_FB(1).orientation(ori).data, [global_min global_max]); % FB-cells D theta
    title('FB B-Cells D theta');
    axis off; colormap(gray); colorbar;

    nexttile;
    imagesc(bPyr2_2_FB(1).orientation(ori).data, [global_min global_max]); % FB-cells D theta+pi
    title('FB B-Cells D theta+pi');
    axis off; colormap(gray); colorbar;

    % Adjust layout size
    set(gcf, 'Position', [100, 100, 1200, 900]); % Adjust figure size
end

% figure; imagesc(params.bPrs.FF.spatial_neighborhood); colorbar;
% figure; imagesc(params.bPrs.FB.spatial_neighborhood); colorbar;