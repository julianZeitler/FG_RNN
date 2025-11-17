close all;
base_path = 'output/feedback/20251002_Overlap_Stimuli';

%% Donut
path = fullfile(base_path, "/donut/all");
params = makeParams("filter_type","donut", "num_scales",3, "G_nu",0.3);

stimulus = ClosedC(700, 400, 320, 192, 64, 1);
runFGSeparation(stimulus, params, true, fullfile(path, "closed C"));

path = fullfile(base_path, "/donut/16");
params = makeParams("filter_type","donut", "num_scales",1, "R0", 16);

stimulus = ClosedC(700, 400, 320, 192, 64, 1);
runFGSeparation(stimulus, params, true, fullfile(path, "closed C"));

path = fullfile(base_path, "/donut/32");
params = makeParams("filter_type","donut", "num_scales",1, "R0", 32);

stimulus = ClosedC(700, 400, 320, 192, 64, 1);
runFGSeparation(stimulus, params, true, fullfile(path, "closed C"));

path = fullfile(base_path, "/donut/64");
params = makeParams("filter_type","donut", "num_scales",1, "R0", 64);

stimulus = ClosedC(700, 400, 320, 192, 64, 1);
runFGSeparation(stimulus, params, true, fullfile(path, "closed C"));

%% Gauss
path = fullfile(base_path, "/gauss/all");
params = makeParams("filter_type","gauss", "num_scales",3, "G_nu",0.3);

stimulus = ClosedC(700, 400, 320, 192, 64, 1);
runFGSeparation(stimulus, params, true, fullfile(path, "closed C"));

path = fullfile(base_path, "/gauss/16");
params = makeParams("filter_type","gauss", "num_scales",1, "R0", 16);

stimulus = ClosedC(700, 400, 320, 192, 64, 1);
runFGSeparation(stimulus, params, true, fullfile(path, "closed C"));

path = fullfile(base_path, "/gauss/32");
params = makeParams("filter_type","gauss", "num_scales",1, "R0", 32);

stimulus = ClosedC(700, 400, 320, 192, 64, 1);
runFGSeparation(stimulus, params, true, fullfile(path, "closed C"));

path = fullfile(base_path, "/gauss/64");
params = makeParams("filter_type","gauss", "num_scales",1, "R0", 64);

stimulus = ClosedC(700, 400, 320, 192, 64, 1);
runFGSeparation(stimulus, params, true, fullfile(path, "closed C"));