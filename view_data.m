groundtruth_path = fullfile('data', 'groundtruth', 'trainfg');
image_path = fullfile('data', 'images');

groundtruth_files = dir(fullfile(groundtruth_path, '*.mat'));

ids = arrayfun(@(x) str2double(erase(x.name, '.mat')), groundtruth_files);

for i = 1:length(ids)
    close all;
    id = ids(i);

    image_filename = fullfile(image_path, 'test', sprintf('%d.jpg', id));
    if ~isfile(image_filename)
        image_filename = fullfile(image_path, 'train', sprintf('%d.jpg', id));
    end
    if ~isfile(image_filename)
       image_filename = fullfile(image_path, 'val', sprintf('%d.jpg', id));
    end
    if ~isfile(image_filename)
        warning(['Image not found for ID: ', num2str(id), ' - skipping image']);
        continue
    end

    disp(['Loaded image: ', image_filename]);

    img = imread(image_filename);

    gt = load(fullfile(groundtruth_path, sprintf('%d.mat', id)));
    gt = gt.gt.groundTruth;
    gt_figure_ground = gt.fg;

    figure; imagesc(img);
    figure; imagesc(gt_figure_ground{1, 1}{1, 1});

    input("\nContinue");
end