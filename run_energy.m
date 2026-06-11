close all;

path = 'output/feedback/20251103_biased_propagation/Long C/';
%% Donut
params = makeParams( ...
    "filter_type","donut", ...
    "num_scales",1, ...
    "G_nu",0.3, ...
    "R0",16, ...
    "iterations",100, ...
    "B_alpha",0.36, ...
    "B_beta",1, ...
    "B_gamma",5, ...
    "B_zeta",0.3, ... 
    "B_lambda",5, ...
    "G_alpha",0.36, ...
    "G_beta",1, ...
    "G_gamma",0, ...
    "G_zeta",0.1, ...
    "G_mu",1);

% stimulus = squareStimulus(256, 256, 128, 1);
% runFGSeparation(stimulus, params, debug=true, log_dir=fullfile(path, "square/1/donut"));
% close all;

stimulus = CStimulus(800, 600, 600, 360, 120, 1);
runFGSeparation(stimulus, params, debug=true, log_dir=fullfile(path, "5/donut"));
close all;

%% Gauss
params = makeParams( ...
    "filter_type","gauss", ...
    "num_scales",1, ...
    "G_nu",0.3, ...
    "R0",16, ...
    "iterations",100, ...
    "B_alpha",0.36, ...
    "B_beta",1, ...
    "B_gamma",5, ...
    "B_zeta",0.3, ...
    "B_lambda",5, ...
    "G_alpha",0.36, ...
    "G_beta",1, ...
    "G_gamma",0, ...
    "G_zeta",0.1, ...
    "G_mu",1);

% stimulus = squareStimulus(256, 256, 128, 1);
% runFGSeparation(stimulus, params, debug=true, log_dir=fullfile(path, "square/1/gauss"));
% close all;

stimulus = CStimulus(800, 600, 600, 360, 120, 1);
runFGSeparation(stimulus, params, debug=true, log_dir=fullfile(path, "5/gauss"));
close all;

% %% Attention
% close all;
% 
% path = 'output/feedback/20251023_attention';
% %% Donut
% params = makeParams("filter_type","donut", "num_scales",3, "G_nu",0.3);
% params.E_scale = 1;
% 
% stimulus = CStimulus(400, 250, 320, 192, 64, 1);
% runFGSeparation(stimulus, params, attention=[0, 1, 0], debug=true, log_dir=fullfile(path, "1/donut"));
% close all;
% 
% runFGSeparation(stimulus, params, attention=[0.05, 0.9, 0.05], debug=true, log_dir=fullfile(path, "2/donut"));
% close all;
% 
% runFGSeparation(stimulus, params, attention=[0.1, 0.8, 0.1], debug=true, log_dir=fullfile(path, "3/donut"));
% close all;
% 
% runFGSeparation(stimulus, params, attention=[0.2, 0.6, 0.2], debug=true, log_dir=fullfile(path, "4/donut"));
% close all;
% 
% runFGSeparation(stimulus, params, attention=[1/3, 1/3, 1/3], debug=true, log_dir=fullfile(path, "5/donut"));
% close all;
% 
% %% Donut
% params = makeParams("filter_type","gauss", "num_scales",3, "G_nu",0.3);
% params.E_scale = 1;
% 
% stimulus = CStimulus(400, 250, 320, 192, 64, 1);
% runFGSeparation(stimulus, params, attention=[0, 1, 0], debug=true, log_dir=fullfile(path, "1/gauss"));
% close all;
% 
% runFGSeparation(stimulus, params, attention=[0.05, 0.9, 0.05], debug=true, log_dir=fullfile(path, "2/gauss"));
% close all;
% 
% runFGSeparation(stimulus, params, attention=[0.1, 0.8, 0.1], debug=true, log_dir=fullfile(path, "3/gauss"));
% close all;
% 
% runFGSeparation(stimulus, params, attention=[0.2, 0.6, 0.2], debug=true, log_dir=fullfile(path, "4/gauss"));
% close all;
% 
% runFGSeparation(stimulus, params, attention=[1/3, 1/3, 1/3], debug=true, log_dir=fullfile(path, "5/gauss"));
% close all;

%%
% log_dir = "output\feedback\20251023_filters_performance\gauss";
% GRF = makeGRF(16,deg2rad([0 22.5 45 67.5 90 112.5 135 157.5])+pi/2, "gauss");
% filter_sum = zeros(size(GRF{1,1}));
% for ori=1:8
%     filter_sum = filter_sum + GRF{1,ori} + GRF{1,ori+8};
%     figure; imagesc(GRF{1,ori}); colorbar; axis image;
%     saveas(gcf, fullfile(log_dir, sprintf('ori_%d.png', ori)));
%     figure; imagesc(GRF{1,ori+8}); colorbar; axis image;
%     saveas(gcf, fullfile(log_dir, sprintf('ori_%d.png', ori+8)));
% end
% figure; imagesc(filter_sum); colorbar; clim([0 max(filter_sum, [], "all")]); axis image;
% saveas(gcf, fullfile(log_dir, 'sum.png'));
% 
% log_dir = "output\feedback\20251023_filters_performance\donut";
% GRF = makeGRF(16,deg2rad([0 22.5 45 67.5 90 112.5 135 157.5])+pi/2, "donut");
% filter_sum = zeros(size(GRF{1,1}));
% for ori=1:8
%     filter_sum = filter_sum + GRF{1,ori} + GRF{1,ori+8};
%     figure; imagesc(GRF{1,ori}); colorbar; axis image;
%     saveas(gcf, fullfile(log_dir, sprintf('ori_%d.png', ori)));
%     figure; imagesc(GRF{1,ori+8}); colorbar; axis image;
%     saveas(gcf, fullfile(log_dir, sprintf('ori_%d.png', ori+8)));
% end
% figure; imagesc(filter_sum); colorbar; clim([0 max(filter_sum, [], "all")]); axis image;
% saveas(gcf, fullfile(log_dir, 'sum.png'));