close all;
params = makeParams();

array = zeros(16, 16);
for g_ori = 1:params.num_ori
    for b_ori = 1:params.num_ori
        array(g_ori,b_ori) = max(-cos(params.oris(g_ori)-params.oris(b_ori)), 0);
        array(g_ori+8,b_ori) = max(-cos(params.oris(g_ori)+pi-params.oris(b_ori)), 0);
        array(g_ori,b_ori+8) = max(-cos(params.oris(g_ori)-params.oris(b_ori)+pi), 0);
        array(g_ori+8,b_ori+8) = max(-cos(params.oris(g_ori)+pi-params.oris(b_ori)+pi), 0);
    end
end
oris = rad2deg([params.oris, params.oris+pi]);
figure; imagesc(array); colorbar;
axis image;
ticks = [1, 3, 5, 7, 9, 11, 13, 15];
set(gca,'YDir','normal', 'XTick', ticks, 'XTickLabel', oris(ticks), 'YTick', ticks, 'YTickLabel', oris(ticks));
fontsize(16, "points")

figure;
tiledlayout(2, 8, 'TileSpacing', 'tight', 'Padding', 'compact');
for ori = 1:length(params.oris)
    nexttile(ori);
    imagesc(params.G.RF{3,ori});
    axis image;
    set(gca,'XTick',[], 'YTick', []);
    nexttile(ori+8);
    imagesc(params.G.RF{3,ori+8});
    axis image;
    set(gca,'XTick',[], 'YTick', []);
end