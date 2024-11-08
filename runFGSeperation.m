function [BOS, edge_map, occ_map, group_map] = runFGSeperation(stimulus, debug, log_dir)
%RUNFGSEPERATION Run the FG Seperation Model on the specified stimulus
%   Inputs:
%       stimulus: 2D numerical array
%       debug: Boolean, if true debug figures will be displayed
%       log_dir: debug figures will be saved to this directory, if
%           specified.s
%
%   Outputs:
%       BOS: mxnx3 Matrix (mxn=model dimensions, 3: HSV-values)
%           H - range [0,1] (encodes [0,2pi]) -> BOS-direction
%           S - normalized strength of BOS-Signal
%           V - 1
%       edge_map: edge detection results
%       occ_map: BOS as image
%       group_map: summed G-Cell activities

debug_flag = false;
save = false;
if nargin >= 2
    if debug==true
        debug_flag = true;
    end

    if nargin >= 3
        save = true;
    end
end

%% Initialize
params = makeParams();
dimensions = size(stimulus);

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

if debug_flag
    % Initialize arrays to store average B activities and BOS-Signals
    avgB_over_iterations = zeros(1, params.iterations);
    BOS_over_iterations = zeros(params.iterations, dimensions(1), dimensions(2), 3);
end

%% Run Model
for idx=1:params.iterations % Loop through all iterations
    [B1, B2] = calculateBCellActivities(E, B1, B2, G, params);
    G = calculateGCellActivities(B1, B2, params);

    if debug_flag
        % Average B1 and B2 over orientations and the entire array (calculate one scalar)
        avgB_over_iterations(idx) = mean(cellfun(@(x) mean(x(:)), {B1})) + mean(cellfun(@(x) mean(x(:)), {B2}));
        % BOS-Signal
        BOS_over_iterations(idx,:,:,:) = getBOS(squeeze(sum(B1,1)), squeeze(sum(B2,1)), params, "unnormalized");
    end
end

%% Output
BOS = getBOS(squeeze(sum(B1,1)), squeeze(sum(B2,1)), params);
edge_map = corfresponse;
occ_map = hsv2rgb(BOS);
group_map = squeeze(sum(G));

%% Debug
if debug_flag
    % 2D graph for average B activity
    figure;
    plot(1:params.iterations, avgB_over_iterations, '-o', 'LineWidth', 2);
    title('Average B Activity Over Iterations');
    xlabel('Iteration');
    ylabel('Average Activity');
    grid on;
    
    if save
        saveas(gcf, fullfile(log_dir, 'average_B_activity.png'));
    end
    
    % Plot BOS over iterations
    BOS_over_iterations(:,:,:,2) = BOS_over_iterations(:,:,:,2)/max(BOS_over_iterations(:,:,:,2), [], "all"); % normalize over all iterations
    figure; tiledlayout(1, params.iterations, 'TileSpacing','compact','Padding','compact');
    for idx=1:params.iterations
        nexttile(idx);
        imagesc(hsv2rgb(squeeze(BOS_over_iterations(idx,:,:,:))));
        title(['Iteration: ', num2str(idx)]);
    end
    set(gcf, "Position", [100,100,320*params.iterations,250]);
    if save
        saveas(gcf, fullfile(log_dir, 'BOS_development.png'));
    end
    
    % Plot BOS on a scale-individual basis
    figure; tiledlayout(2, 5, 'TileSpacing', 'compact', 'Padding', 'compact');
    for k=1:params.num_scales
        nexttile(k);
        imagesc(hsv2rgb(getBOS(squeeze(B1(k,:,:,:)), squeeze(B2(k,:,:,:)), params)));
        title(['G-Cell RF Radius = ', num2str(params.scale_step^k)]);
    end
    set(gcf, 'Position', [50, 50, 1800, 600]);
    if save
        saveas(gcf, fullfile(log_dir, 'BOS_scale-individual.png'));
    end
    
    % Plot G-activity on a scale-individual basis
    figure; tiledlayout(2, 5, 'TileSpacing', 'compact', 'Padding', 'tight');
    for k=1:params.num_scales
        nexttile(k);
        imagesc(squeeze(G(k,:,:))); colorbar;
        title(['G-Cell RF Radius = ', num2str(params.scale_step^k)]);
    end
    set(gcf, 'Position', [50, 50, 1900, 600]);
    if save
        saveas(gcf, fullfile(log_dir, 'G_scale-individual.png'));
    end
    
    % Summed G-Cell activity
    figure; imagesc(squeeze(sum(G))); colorbar;
    if save
        saveas(gcf, fullfile(log_dir, 'summed_G.png'));
    end
    
    % Collapsed Scales BOS-Signal
    figure; imagesc(hsv2rgb(getBOS(squeeze(sum(B1,1)), squeeze(sum(B2,1)), params)));
    if save
        saveas(gcf, fullfile(log_dir, 'BOS.png'));
    end
    
    % Stimulus
    figure; imagesc(stimulus); colormap(gray); colorbar;
    if save
        saveas(gcf, fullfile(log_dir, 'stimulus.png'));
    end
    
    % Edge Detection
    figure; imagesc(corfresponse); axis image; colormap(gray); colorbar;
    if save
        saveas(gcf, fullfile(log_dir, 'edge_detection.png'));
    end
    
    % B-Cell activity and BOS-Signal split for orientation
    figure; 
    tiledlayout(3, 8, 'TileSpacing', 'tight', 'Padding', 'tight');
    % Plot B and BOS for each orientation
    % Loop over orientations (1 to 8)
    for ori = 1:8
        % Plot B1 cells for each orientation in the first row
        nexttile(ori); % First row
        imagesc(squeeze(sum(B1(:,ori,:,:))));
        title(['B1 \theta=', num2str(ori)]);
        colorbar;
    
        % Plot B2 cells for each orientation in the second row
        nexttile(ori + 8); % Second row
        imagesc(squeeze(sum(B2(:,ori,:,:))));
        title(['B2 \theta=', num2str(ori)]);
        colorbar;
    
        % Plot BOS signal for each orientation in the third row
        nexttile(ori + 16); % Third row
        BOS_signal = squeeze(sum(B1(:,ori,:,:)))-squeeze(sum(B2(:,ori,:,:)));
        imagesc(BOS_signal);
        title(['BOS \theta=', num2str(ori)]);
        colorbar;
    end
    set(gcf, 'Position', [50, 50, 1800, 500]);
    if save
        saveas(gcf, fullfile(log_dir, 'B_BOS.png'));
    end
end
end

