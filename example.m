close all;

%% Natural Stimuli
figure;
tiledlayout(3, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

stimulus = rgb2gray(imread("images/42049.jpg"));
[BOS, edge_map, occ_map, group_map] = runFGSeparation(stimulus);
group_map = group_map./max(group_map, [], "all");
ax1 = nexttile;
imagesc(stimulus); axis image; set(gca,'XTick',[], 'YTick', []); colormap(ax1, "gray");
ax2 = nexttile;
imagesc(edge_map); axis image; set(gca,'XTick',[], 'YTick', []); colormap(ax2, "gray");
nexttile;
imagesc(occ_map); axis image; set(gca,'XTick',[], 'YTick', []);
ax3 = nexttile;
imagesc(group_map); colorbar; axis image; set(gca,'XTick',[], 'YTick', []); colormap(ax3, "parula");

stimulus = rgb2gray(imread("images/12074.jpg"));
[BOS, edge_map, occ_map, group_map] = runFGSeparation(stimulus);
ax1 = nexttile;
imagesc(stimulus); axis image; set(gca,'XTick',[], 'YTick', []); colormap(ax1, "gray");
ax2 = nexttile;
imagesc(edge_map); axis image; set(gca,'XTick',[], 'YTick', []); colormap(ax2, "gray");
nexttile;
imagesc(occ_map); axis image; set(gca,'XTick',[], 'YTick', []);
ax3 = nexttile;
imagesc(group_map); axis image; set(gca,'XTick',[], 'YTick', []); colormap(ax3, "parula");

stimulus = rgb2gray(imread("images/156079.jpg"));
[BOS, edge_map, occ_map, group_map] = runFGSeparation(stimulus);
ax1 = nexttile;
imagesc(stimulus); axis image; set(gca,'XTick',[], 'YTick', []); colormap(ax1, "gray");
ax2 = nexttile;
imagesc(edge_map); axis image; set(gca,'XTick',[], 'YTick', []); colormap(ax2, "gray");
nexttile;
imagesc(occ_map); axis image; set(gca,'XTick',[], 'YTick', []);
ax3 = nexttile;
imagesc(group_map); axis image; set(gca,'XTick',[], 'YTick', []); colormap(ax3, "parula");

%% Non-Natural Stimuli
dimensions = [200 200];
figure;
tiledlayout(4, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

stimulus = squareStimulus(dimensions(1), dimensions(2), 50, 1);
[BOS, edge_map, occ_map, group_map] = runFGSeparation(stimulus);
group_map = group_map./max(group_map, [], "all");
ax1 = nexttile;
imagesc(stimulus); axis image; set(gca,'XTick',[], 'YTick', []); colormap(ax1, "gray");
ax2 = nexttile;
imagesc(edge_map); axis image; set(gca,'XTick',[], 'YTick', []); colormap(ax2, "gray");
nexttile;
imagesc(occ_map); axis image; set(gca,'XTick',[], 'YTick', []);
ax3 = nexttile;
imagesc(group_map); colorbar; axis image; set(gca,'XTick',[], 'YTick', []); colormap(ax3, "parula");

stimulus = changingBackground(dimensions(1), dimensions(2), 50);
[BOS, edge_map, occ_map, group_map] = runFGSeparation(stimulus);
ax1 = nexttile;
imagesc(stimulus); axis image; set(gca,'XTick',[], 'YTick', []); colormap(ax1, "gray");
ax2 = nexttile;
imagesc(edge_map); axis image; set(gca,'XTick',[], 'YTick', []); colormap(ax2, "gray");
nexttile;
imagesc(occ_map); axis image; set(gca,'XTick',[], 'YTick', []);
ax3 = nexttile;
imagesc(group_map); axis image; set(gca,'XTick',[], 'YTick', []); colormap(ax3, "parula");

stimulus = vaseStimulus(dimensions(1), dimensions(2), 5, 50);
[BOS, edge_map, occ_map, group_map] = runFGSeparation(stimulus);
ax1 = nexttile;
imagesc(stimulus); axis image; set(gca,'XTick',[], 'YTick', []); colormap(ax1, "gray");
ax2 = nexttile;
imagesc(edge_map); axis image; set(gca,'XTick',[], 'YTick', []); colormap(ax2, "gray");
nexttile;
imagesc(occ_map); axis image; set(gca,'XTick',[], 'YTick', []);
ax3 = nexttile;
imagesc(group_map); axis image; set(gca,'XTick',[], 'YTick', []); colormap(ax3, "parula");

stimulus = CStimulus(dimensions(1), 50, 10, 1);
[BOS, edge_map, occ_map, group_map] = runFGSeparation(stimulus);
ax1 = nexttile;
imagesc(stimulus); axis image; set(gca,'XTick',[], 'YTick', []); colormap(ax1, "gray");
ax2 = nexttile;
imagesc(edge_map); axis image; set(gca,'XTick',[], 'YTick', []); colormap(ax2, "gray");
nexttile;
imagesc(occ_map); axis image; set(gca,'XTick',[], 'YTick', []);
ax3 = nexttile;
imagesc(group_map); axis image; set(gca,'XTick',[], 'YTick', []); colormap(ax3, "parula");
