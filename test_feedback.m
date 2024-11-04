close all; clear;
params = makeParams();
stimulus = rgb2gray(imread("images\42049.jpg"));
dimensions = size(stimulus);
% dimensions = [200 200];
% stimulus = squareStimulus(dimensions(1), dimensions(2), 50, 1);
% stimulus = verticalContrastStimulus(dimensions(1), dimensions(2));
% stimulus = changingBackground(dimensions(1), dimensions(2), 50);
% stimulus = verticalBarStimulus(dimensions(1), dimensions(2), 50);
% stimulus = vaseStimulus(200, 200, 5, 50);
% stimulus = CStimulus(200,50,10,1);

% save = true;
save = false;

% Initialize activities to zero

B1 = zeros(params.num_scales, params.num_ori, dimensions(1), dimensions(2));
B2 = zeros(params.num_scales, params.num_ori, dimensions(1), dimensions(2));
E = zeros(params.num_ori, dimensions(1), dimensions(2));
G = zeros(params.num_scales, dimensions(1), dimensions(2));

% corfresponse contains contrasts, oriensMatrix contains the orientation of
% contrasts at each spatial location
[~,~,corfresponse,oriensMatrix] = CORFContourDetection(stimulus,2.2,4,1.8);
oriensMatrix = mod(round(oriensMatrix/(2*pi)*16)+4,16)+1; % convert orientations to indices

for ori=1:params.num_ori
    % Combine opposite contrast polarities at each orientation
    E(ori, :, :) = corfresponse.*(oriensMatrix==ori) + corfresponse.*(oriensMatrix==ori+8);
end

% Initialize arrays to store average B activities and BOS-Signals
avgB_over_iterations = zeros(1, params.iterations);
BOS_over_iterations = zeros(params.iterations, dimensions(1), dimensions(2), 3);

for idx=1:params.iterations % Loop through all iterations
    [B1, B2] = calculateBCellActivities(E, B1, B2, G, params);
    
    G = calculateGCellActivities(B1, B2, params);

    % Average B1 and B2 over orientations and the entire array (calculate one scalar)
    avgB_over_iterations(idx) = mean(cellfun(@(x) mean(x(:)), {B1})) + mean(cellfun(@(x) mean(x(:)), {B2}));
    % BOS-Signal
    BOS_over_iterations(idx,:,:,:) = getBOS(squeeze(sum(B1,1)), squeeze(sum(B2,1)), params, "unnormalized");
end

%% Visualize
% 2D graph for average B activity
figure;
plot(1:params.iterations, avgB_over_iterations, '-o', 'LineWidth', 2);
title('Average B Activity Over Iterations');
xlabel('Iteration');
ylabel('Average Activity');
grid on;

if save
    saveas(gcf, 'output/feedback/average_B_activity.png');
end

% Plot BOS over iterations
BOS_over_iterations(:,:,:,2) = BOS_over_iterations(:,:,:,2)/max(BOS_over_iterations(:,:,:,2), [], "all"); % normalize over all iterations
figure; tiledlayout(1, params.iterations, 'TileSpacing','tight','Padding','tight');
for idx=1:params.iterations
    nexttile(idx);
    imagesc(hsv2rgb(squeeze(BOS_over_iterations(idx,:,:,:))));
    title(['Iteration: ', num2str(idx)]);
    axis off;
end
set(gcf, "Position", [100,100,3000,250]);
if save
    saveas(gcf, 'output/feedback/BOS_development.png');
end

% Miscellaneous
figure; imagesc(squeeze(sum(G))); colorbar; axis off;
if save
    saveas(gcf, 'output/feedback/summed_G.png');
end
figure; imagesc(hsv2rgb(getBOS(squeeze(sum(B1,1)), squeeze(sum(B2,1)), params)));
if save
    saveas(gcf, 'output/feedback/BOS.png');
end
% figure; imagesc(stimulus); colormap(gray); axis off; colorbar;
% figure; imagesc(corfresponse); axis image; axis off; colormap(gray); colorbar;

figure; 
tiledlayout(3, 8, 'TileSpacing', 'tight', 'Padding', 'tight');

% Plot B and BOS for each iteration
% Loop over orientations (1 to 8)
for ori = 1:8
    % Plot B1 cells for each orientation in the first row
    nexttile(ori); % First row
    imagesc(squeeze(sum(B1(:,ori,:,:))));
    title(['B1 \theta=', num2str(ori)]);
    colorbar; axis off;

    % Plot B2 cells for each orientation in the second row
    nexttile(ori + 8); % Second row
    imagesc(squeeze(sum(B2(:,ori,:,:))));
    title(['B2 \theta=', num2str(ori)]);
    colorbar; axis off;

    % Plot BOS signal for each orientation in the third row
    nexttile(ori + 16); % Third row
    BOS_signal = squeeze(sum(B1(:,ori,:,:)))-squeeze(sum(B2(:,ori,:,:)));
    imagesc(BOS_signal);
    title(['BOS \theta=', num2str(ori)]);
    colorbar; axis off;
end

set(gcf, 'Position', [50, 50, 1800, 500]);
if save
    saveas(gcf, 'output/feedback/B_BOS.png');
end