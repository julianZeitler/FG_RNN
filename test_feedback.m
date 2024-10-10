close all; clear;
params = makeParams();
stimulus = rgb2gray(imread("images\12074.jpg"));
dimensions = size(stimulus);
% dimensions = [200 200];
% stimulus = squareStimulus(dimensions(1), dimensions(2), 50, 1);
% stimulus = verticalContrastStimulus(dimensions(1), dimensions(2));
% stimulus = changingBackground(dimensions(1), dimensions(2), 50);
% stimulus = verticalBarStimulus(dimensions(1), dimensions(2), 50);
% stimulus = vaseStimulus(200, 200, 5, 50);

% Initialize activities to zero

B1 = zeros(params.numOri, dimensions(1), dimensions(2));
B2 = zeros(params.numOri, dimensions(1), dimensions(2));
E = zeros(params.numOri, dimensions(1), dimensions(2));
G = zeros(dimensions(1), dimensions(2));

% corfresponse contains contrasts, oriensMatrix contains the orientation of
% contrasts at each spatial location
[~,~,corfresponse,oriensMatrix] = CORFContourDetection(stimulus,2.2,4,1.8);
oriensMatrix = mod(round(oriensMatrix/(2*pi)*16)+4,16)+1; % convert orientations to indices

for ori=1:params.numOri
    % Combine opposite contrast polarities at each orientation
    E(ori, :, :) = corfresponse.*(oriensMatrix==ori) + corfresponse.*(oriensMatrix==ori+8);
end

% Original figure for G-cells, B1 (orientation 5), and B2 (orientation 5)
% figure; tiledlayout(3, params.iterations, 'TileSpacing', 'tight', 'Padding', 'tight'); 

% Separate figures for 2D plot of averaged B1 and B2 activities
figureB1 = figure; hold on;
figureB2 = figure; hold on;

% Initialize arrays to store average B1 and B2 activities
avgB1_over_iterations = zeros(1, params.iterations);
avgB2_over_iterations = zeros(1, params.iterations);

for idx=1:params.iterations % Loop through all iterations
    [B1, B2] = calculateBCellActivities(E, B1, B2, G, params);
    
    G = calculateGCellActivities(B1, B2, params);

    % Average B1 and B2 over orientations and the entire array (calculate one scalar)
    avgB1_over_iterations(idx) = mean(cellfun(@(x) mean(x(:)), {B1}));
    avgB2_over_iterations(idx) = mean(cellfun(@(x) mean(x(:)), {B2}));

end

% Additional Plot: 2D graph for average B1 activity
figure(figureB1);
plot(1:params.iterations, avgB1_over_iterations, '-o', 'LineWidth', 2);
title('Average B1 Activity Over Iterations');
xlabel('Iteration');
ylabel('Average Activity');
grid on;

% Additional Plot: 2D graph for average B2 activity
figure(figureB2);
plot(1:params.iterations, avgB2_over_iterations, '-o', 'LineWidth', 2);
title('Average B2 Activity Over Iterations');
xlabel('Iteration');
ylabel('Average Activity');
grid on;

figure; imagesc(G); colorbar; axis off;

%% Visualize
figure; imagesc(stimulus); colormap(gray); axis off; colorbar;
figure; imagesc(corfresponse); axis image; axis off; colormap(gray); colorbar;

figure; 
tiledlayout(3, 8, 'TileSpacing', 'tight', 'Padding', 'tight');

% Loop over orientations (1 to 8)
for ori = 1:8
    % Plot B1 cells for each orientation in the first row
    nexttile(ori); % First row
    imagesc(squeeze(B1(ori,:,:)));
    title(['B1 \theta=', num2str(ori)]);
    colorbar; axis off;
    
    % Plot B2 cells for each orientation in the second row
    nexttile(ori + 8); % Second row
    imagesc(squeeze(B2(ori,:,:)));
    title(['B2 \theta=', num2str(ori)]);
    colorbar; axis off;
    
    % Plot BOS signal for each orientation in the third row
    nexttile(ori + 16); % Third row
    BOS_signal = squeeze(B1(ori,:,:)) - squeeze(B2(ori,:,:));
    imagesc(BOS_signal);
    title(['BOS \theta=', num2str(ori)]);
    colorbar; axis off;
end

set(gcf, 'Position', [50, 50, 1800, 500]);