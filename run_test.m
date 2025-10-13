close all;

base_path = 'output/feedback/20251001_GRFs';
%% Donut
path = fullfile(base_path, "/donut/3_scales_competition");
params = makeParams("filter_type","donut", "num_scales",3, "G_nu",0.3);

stimulus = CStimulus(400, 250, 320, 192, 64, 1);
runFGSeparation(stimulus, params, true, fullfile(path, "C"));
close all;

stimulus = overlappingRectangles(256, 256, 192, 128, 64, 64);
runFGSeparation(stimulus, params, true, fullfile(path, "overlap"));
close all;

stimulus = squareStimulus(256, 256, 64, 1);
runFGSeparation(stimulus, params, true, fullfile(path, "square"));
close all;

stimulus = rgb2gray(imread("images/42049.jpg"));
runFGSeparation(stimulus, params, true, fullfile(path, "eagle"));
close all;

stimulus = rgb2gray(imread("images/12074.jpg"));
runFGSeparation(stimulus, params, true, fullfile(path, "mushroom"));
close all;


path = fullfile(base_path, "/donut/1_scale/32");
params = makeParams("filter_type","donut", "num_scales",1, "R0", 32);

stimulus = rgb2gray(imread("images/42049.jpg"));
runFGSeparation(stimulus, params, true, fullfile(path, "eagle"));
close all;

stimulus = rgb2gray(imread("images/12074.jpg"));
runFGSeparation(stimulus, params, true, fullfile(path, "mushroom"));
close all;

stimulus = CStimulus(400, 250, 320, 192, 64, 1);
runFGSeparation(stimulus, params, true, fullfile(path, "C"));
close all;

stimulus = overlappingRectangles(256, 256, 192, 128, 64, 64);
runFGSeparation(stimulus, params, true, fullfile(path, "overlap"));
close all;

stimulus = overlappingRectangles(320, 320, 192, 192, 64, 128);
runFGSeparation(stimulus, params, true, fullfile(path, "overlap C"));
close all;

stimulus = overlappingRectangles(256, 256, 192, 192, 64, 64);
runFGSeparation(stimulus, params, true, fullfile(path, "overlap big C"));
close all;

stimulus = squareStimulus(256, 256, 64, 1);
runFGSeparation(stimulus, params, true, fullfile(path, "square"));
close all;

% Additional plots for overlap stimulus
path = fullfile(base_path, "/donut/1_scale/16");
params = makeParams("filter_type","donut", "num_scales",1, "R0",16);

stimulus = overlappingRectangles(256, 256, 192, 128, 64, 64);
runFGSeparation(stimulus, params, true, fullfile(path, "overlap"));
close all;

stimulus = overlappingRectangles(320, 320, 192, 192, 64, 128);
runFGSeparation(stimulus, params, true, fullfile(path, "overlap C"));
close all;

stimulus = overlappingRectangles(256, 256, 192, 192, 64, 64);
runFGSeparation(stimulus, params, true, fullfile(path, "overlap big C"));
close all;

path = fullfile(base_path, "/donut/1_scale/64");
params = makeParams("filter_type","donut", "num_scales",1, "R0",64);

stimulus = overlappingRectangles(256, 256, 192, 128, 64, 64);
runFGSeparation(stimulus, params, true, fullfile(path, "overlap"));
close all;

stimulus = overlappingRectangles(320, 320, 192, 192, 64, 128);
runFGSeparation(stimulus, params, true, fullfile(path, "overlap C"));
close all;

stimulus = overlappingRectangles(256, 256, 192, 192, 64, 64);
runFGSeparation(stimulus, params, true, fullfile(path, "overlap big C"));
close all;


%% Gauss
path = fullfile(base_path, "/gauss/3_scales_competition_fixed_sigma");
params = makeParams("filter_type","gauss", "num_scales",3,"G_nu",0.3);

stimulus = CStimulus(400, 250, 320, 192, 64, 1);
runFGSeparation(stimulus, params, true, fullfile(path, "C"));
close all;

stimulus = overlappingRectangles(256, 256, 192, 128, 64, 64);
runFGSeparation(stimulus, params, true, fullfile(path, "overlap"));
close all;

stimulus = squareStimulus(256, 256, 64, 1);
runFGSeparation(stimulus, params, true, fullfile(path, "square"));
close all;

stimulus = rgb2gray(imread("images/42049.jpg"));
runFGSeparation(stimulus, params, true, fullfile(path, "eagle"));
close all;

stimulus = rgb2gray(imread("images/12074.jpg"));
runFGSeparation(stimulus, params, true, fullfile(path, "mushroom"));
close all;


path = fullfile(base_path, "/gauss/1_scale_fixed_gauss/32");
params = makeParams("filter_type","gauss", "num_scales",1, "R0",32);

stimulus = rgb2gray(imread("images/42049.jpg"));
runFGSeparation(stimulus, params, true, fullfile(path, "eagle"));
close all;

stimulus = rgb2gray(imread("images/12074.jpg"));
runFGSeparation(stimulus, params, true, fullfile(path, "mushroom"));
close all;

stimulus = CStimulus(400, 250, 320, 192, 64, 1);
runFGSeparation(stimulus, params, true, fullfile(path, "C"));
close all;

stimulus = overlappingRectangles(256, 256, 192, 128, 64, 64);
runFGSeparation(stimulus, params, true, fullfile(path, "overlap"));
close all;

stimulus = overlappingRectangles(320, 320, 192, 192, 64, 128);
runFGSeparation(stimulus, params, true, fullfile(path, "overlap C"));
close all;

stimulus = overlappingRectangles(256, 256, 192, 192, 64, 64);
runFGSeparation(stimulus, params, true, fullfile(path, "overlap big C"));
close all;

stimulus = squareStimulus(256, 256, 64, 1);
runFGSeparation(stimulus, params, true, fullfile(path, "square"));
close all;

% Additional plots for overlap stimulus
path = fullfile(base_path, "/gauss/1_scale_fixed_gauss/16");
params = makeParams("filter_type","gauss", "num_scales",1, "R0",16);

stimulus = overlappingRectangles(256, 256, 192, 128, 64, 64);
runFGSeparation(stimulus, params, true, fullfile(path, "overlap"));
close all;

stimulus = overlappingRectangles(320, 320, 192, 192, 64, 128);
runFGSeparation(stimulus, params, true, fullfile(path, "overlap C"));
close all;

stimulus = overlappingRectangles(256, 256, 192, 192, 64, 64);
runFGSeparation(stimulus, params, true, fullfile(path, "overlap big C"));
close all;

path = fullfile(base_path, "/gauss/1_scale_fixed_gauss/64");
params = makeParams("filter_type","gauss", "num_scales",1, "R0",64);

stimulus = overlappingRectangles(256, 256, 192, 128, 64, 64);
runFGSeparation(stimulus, params, true, fullfile(path, "overlap"));
close all;

stimulus = overlappingRectangles(320, 320, 192, 192, 64, 128);
runFGSeparation(stimulus, params, true, fullfile(path, "overlap C"));
close all;

stimulus = overlappingRectangles(256, 256, 192, 192, 64, 64);
runFGSeparation(stimulus, params, true, fullfile(path, "overlap big C"));
close all;



params = makeParams("filter_type","gauss", "num_scales",3, "G_nu",0.3);
runFGSeparation(stimulus, params, true, fullfile("output/feedback/20251008_scale_compensation/gauss/attention_1/C"));
close all;

params = makeParams("filter_type","donut", "num_scales",3, "G_nu",0.3);
runFGSeparation(stimulus, params, true, fullfile("output/feedback/20251008_scale_compensation/donut/attention_1/C"));
close all;