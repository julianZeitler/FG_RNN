close all;

savePath = 'G:\Meine Ablage\Studium\Master\Project Modulatory Feedback\FG_RNN\output\feedback\20250918_back_to_2D\Line';
if ~exist(savePath, 'dir')
    mkdir(savePath);
end

params.iterations = 50;
params.oris = deg2rad([0, 90]);
params.num_ori = length(params.oris);

params.B.alpha = 0.36;
params.B.beta = 1;
params.B.gamma = 7;
params.B.zeta = 0.7;
params.B.lambda = 5;

params.G.alpha = 0.25;  
params.G.beta = 1;
params.G.gamma = 0;
params.G.zeta = 0.9;
params.G.mu = 1; % inhibition from opposing B cells
params.G.inhibition_filter = gaussianFilter2D(5, 5, 1, 1);
params.G.RF = makeGRF(1, params.oris+pi/2);

[G_loc, BV_loc, BH_loc, E_loc] = parseStringMap(getStringMap(5));

B_loc = zeros(params.num_ori, size(BH_loc, 1), size(BH_loc, 2));
B_loc(1, :, :) = BH_loc;
B_loc(2, :, :) = BV_loc;
E = zeros(params.num_ori, size(BH_loc, 1), size(BH_loc, 2));
E(1, :, :) = E_loc.*squeeze(B_loc(1, :, :));
E(2, :, :) = E_loc.*squeeze(B_loc(2, :, :));

G = zeros(size(G_loc));
B1 = zeros(size(B_loc));
B2 = zeros(size(B_loc));

B1_prev = B1;
B2_prev = B2;
B1_history = zeros([size(B1_prev), params.iterations+1]);
B2_history = zeros([size(B2_prev), params.iterations+1]);
G_history = zeros([size(G), params.iterations+1]);

B1_history(:,:,:,1) = B1_prev;
B2_history(:,:,:,1) = B2_prev;
G_history(:,:,1) = G;

BOS_over_iterations = zeros([size(G), 3, params.iterations]);

for i = 1:params.iterations
    for ori = 1:params.num_ori
        FB1 = imfilter(G, params.G.RF{ori+params.num_ori});
        B1(ori, :, :) = params.B.beta*squeeze(E(ori,:,:)).*(1 + params.B.lambda*FB1)./(...
            params.B.alpha + ...
            params.B.gamma*squeeze(B2_prev(ori,:,:)) + ...
            params.B.zeta*squeeze(E(ori,:,:)).*(1 + params.B.lambda*FB1));


        FB2 = imfilter(G, params.G.RF{ori});
        B2(ori, :, :) = params.B.beta*squeeze(E(ori,:,:)).*(1 + params.B.lambda*FB2)./(...
            params.B.alpha + ...
            params.B.gamma*squeeze(B1_prev(ori,:,:)) + ...
            params.B.zeta*squeeze(E(ori,:,:)).*(1 + params.B.lambda*FB2));
    end
    G_in = zeros(size(G));
    G_in_inh = zeros(size(G));
    for ori = 1:params.num_ori
        G_in = G_in + imfilter(squeeze(B1(ori,:,:)), params.G.RF{ori});
        G_in = G_in + imfilter(squeeze(B2(ori,:,:)), params.G.RF{ori+params.num_ori});

        G_in_inh = G_in_inh + imfilter(squeeze(B1(ori,:,:)), params.G.RF{ori+params.num_ori});
        G_in_inh = G_in_inh + imfilter(squeeze(B2(ori,:,:)), params.G.RF{ori});
    end
    G = params.G.beta*G_in./( ...
        params.G.alpha + ...
        params.G.gamma*imfilter(G, params.G.inhibition_filter) + ...
        params.G.zeta*G_in + ...
        params.G.mu*G_in_inh);

    B1_history(:,:,:,i+1) = B1;
    B2_history(:,:,:,i+1) = B2;
    G_history(:,:,i+1)   = G;

    B1_prev = B1;
    B2_prev = B2;

    BOS_over_iterations(:,:,:,i) = getBOS(B1, B2, params, "unnormalized");
end


% figure; imagesc(E); axis image; colorbar;
% figure; imagesc(BV1_loc); axis image;
% figure; imagesc(BH1_loc); axis image;
% figure; imagesc(G_loc); axis image;
figure; imagesc(hsv2rgb(getBOS(B1, B2, params))); axis image;
saveas(gcf, fullfile(savePath, 'BOS.png'));


horizontal_locations = [2,2; 2,4; 2,6; 2,8; 2,16];
vertical_locations = [1,1;];
g_locations = [1,2; 1,4; 1,8; 1,16];

numVLocations = size(vertical_locations, 1);
numHLocations = size(horizontal_locations, 1);
numGLocations = size(g_locations, 1);
maxPerRow = 5;

numVRows = ceil(numVLocations / maxPerRow);
numHRows = ceil(numHLocations / maxPerRow);
numGRows = ceil(numGLocations / maxPerRow);
iterations_axis = 1:params.iterations+1;

%% --- Vertical B Cells ---
figure('Position', [100, 100, 1200, 600]);
for idx = 1:numVLocations
    row = vertical_locations(idx, 1);
    col = vertical_locations(idx, 2);

    BV1_trace = squeeze(B1_history(2, row, col, :));
    BV2_trace = squeeze(B2_history(2, row, col, :));

    subplot(numVRows, min(maxPerRow, numVLocations), idx);
    plot(iterations_axis, BV1_trace, 'b-', 'LineWidth', 1.5); hold on;
    plot(iterations_axis, BV2_trace, 'r-', 'LineWidth', 1.5);
    title(sprintf('Vertical (%d, %d)', row, col));
    xlabel('Iter'); ylabel('Act.');
    legend('BV1','BV2','Location','southeast');
    grid on;
end
sgtitle('Vertical B Cell Activity Over Iterations');
set(gcf, "Position", [100, 100, 1200, 300*numVRows]);
saveas(gcf, fullfile(savePath, 'Vertical B Cells.png'));

%% --- Horizontal B Cells ---
figure;
for idx = 1:numHLocations
    row = horizontal_locations(idx, 1);
    col = horizontal_locations(idx, 2);

    BH1_trace = squeeze(B1_history(1, row, col, :));
    BH2_trace = squeeze(B2_history(1, row, col, :));

    subplot(numHRows, min(maxPerRow, numHLocations), idx);
    plot(iterations_axis, BH1_trace, 'b-', 'LineWidth', 1.5); hold on;
    plot(iterations_axis, BH2_trace, 'r-', 'LineWidth', 1.5);
    title(sprintf('Horizontal (%d, %d)', row, col));
    xlabel('Iter'); ylabel('Act.');
    legend('BH1','BH2','Location','southeast');
    grid on;
end
sgtitle('Horizontal B Cell Activity Over Iterations');
set(gcf, "Position", [100, 100, 1200, 300*numHRows]);
saveas(gcf, fullfile(savePath, 'Horizontal B Cells.png'));

%% --- G Cells ---
figure;
for idx = 1:numGLocations
    row = g_locations(idx, 1);
    col = g_locations(idx, 2);

    G_trace = squeeze(G_history(row, col, :));

    subplot(numGRows, min(maxPerRow, numGLocations), idx);
    plot(iterations_axis, G_trace, 'k-', 'LineWidth', 1.5);
    title(sprintf('G Cell (%d, %d)', row, col));
    xlabel('Iter'); ylabel('Act.');
    legend('G','Location','southeast');
    grid on;
end
sgtitle('G Cell Activity Over Iterations');
set(gcf, "Position", [100, 100, 1200, 300*numGRows]);
saveas(gcf, fullfile(savePath, 'G Cells.png'));

%% BOS over iterations
BOS_over_iterations(:,:,:,2) = BOS_over_iterations(:,:,:,2)/max(BOS_over_iterations(:,:,:,2), [], "all");
figure('Position', [100, 100, 1500, 300]);
selected_iterations = round(linspace(1, params.iterations, 5));
for idx = 1:5
    iter = selected_iterations(idx);
    subplot(1, 5, idx);
    imagesc(hsv2rgb(BOS_over_iterations(:,:,:,iter))); 
    axis image;
    title(sprintf('Iteration %d', iter));
    axis off;
end
sgtitle('BOS Evolution Over Iterations');
saveas(gcf, fullfile(savePath, 'BOS_over_iterations.png'));