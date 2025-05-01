close all;
path = 'output/feedback/20250430_Scale_Competition';

stimulus = rgb2gray(imread("images/42049.jpg"));
runFGSeparation(stimulus, true, fullfile(path, "eagle"));
close all;

stimulus = rgb2gray(imread("images/12074.jpg"));
runFGSeparation(stimulus, true, fullfile(path, "mushroom"));
close all;

dimensions = [200 200];

stimulus = squareStimulus(dimensions(1), dimensions(2), 50, 1);
runFGSeparation(stimulus, true, fullfile(path, "square"));
close all;

stimulus = changingBackground(dimensions(1), dimensions(2), 50);
runFGSeparation(stimulus, true, fullfile(path, "gradient"));
close all;

stimulus = vaseStimulus(dimensions(1), dimensions(2), 5, 50);
runFGSeparation(stimulus, true, fullfile(path, "vase"));
close all;

stimulus = CStimulus(dimensions(1), 50, 10, 1);
runFGSeparation(stimulus, true, fullfile(path, "c"));
close all;

