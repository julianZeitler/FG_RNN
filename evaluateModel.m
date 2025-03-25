groundtruth_path = fullfile('data', 'groundtruth', 'trainfg');
image_path = fullfile('data', 'images');

groundtruth_files = dir(fullfile(groundtruth_path, '*.mat'));

% Extract IDs from filenames
ids = arrayfun(@(x) str2double(erase(x.name, '.mat')), groundtruth_files);

AccF = zeros(length(ids), 1);
diffS = zeros(length(ids), 1);

for i = 1:length(ids)
    id = ids(i); 

    % Construct the filename for the image id.jpg. Look in all image
    % directories
    image_filename = fullfile(image_path, 'test', sprintf('%d.jpg', id));
    if ~isfile(image_filename)
        image_filename = fullfile(image_path, 'train', sprintf('%d.jpg', id));
    end
    if ~isfile(image_filename)
       image_filename = fullfile(image_path, 'val', sprintf('%d.jpg', id));
    end
    if ~isfile(image_filename)
        warning(['Image not found for ID: ', num2str(id), ' - skipping image']);
        AccF(i,1) = NaN;
        diffS(i,1) = NaN;
        continue
    end

    % Prediction
    img = rgb2gray(imread(image_filename));
    disp(['Loaded image: ', image_filename]);

    [BOS, edge_map, occ_map, group_map] = runFGSeparation(img);
    BOS_strength = squeeze(BOS(:,:,2));
    BOS_orientation = 2*pi*squeeze(BOS(:,:,1)); % Orientations are between 0 and 1

    % Ground Truth
    gt = load(fullfile(groundtruth_path, sprintf('%d.mat', id)));
    gt = gt.gt.groundTruth;
    gt_figure_ground = gt.fg;
    n_figure_ground_assignments = size(gt_figure_ground,2);

    gt_foreground = zeros(size(occ_map,1), size(occ_map,2), n_figure_ground_assignments);
    gt_background = gt_foreground;
     
    for k=1:n_figure_ground_assignments
        gt_foreground(:,:,k) = gt_figure_ground{1,k}{1,1}==1;
        gt_background(:,:,k) = gt_figure_ground{1,k}{1,1}==-1;
    end

    % 
    gt_foreground=any(gt_foreground,3);
    gt_background=any(gt_background,3);
    
    % expand and combine -- get 1 pixel width
    gt_foreground_collapsed = gt_foreground;
    gt_foreground_collapsed(bwdist(gt_foreground) <= 1) = 1;
    gt_background_collapsed = gt_background;
    gt_background_collapsed(bwdist(gt_background) <= 1) = 1;
    gt_combined = and(gt_foreground_collapsed, gt_background_collapsed);
    
    % combine GT and GTMasks
    gt_combined_edges = zeros(size(gt_foreground));
    gt_combined_edges(gt_foreground) = 1;
    gt_combined_edges(gt_background) = -1;
    gt_combined_edges(~gt_combined) = 0;
    
    % convert groundtruth_combined_edges to orientation code (1 pixel)
    [~,gt_ori] = computeEdge3(gt_combined_edges, gt_combined_edges>0);
    gt_ori = wrapTo2Pi(gt_ori-pi/2);
    
    [M1,~] = correspondPixels(double(BOS_strength>0.01), double(gt_ori>0), 0.0075);
    iResP = find(M1);
    igtP = M1(find(M1));
    iOC = [BOS_orientation(iResP), gt_ori(igtP)];
    diffC = abs(iOC(:,1) - iOC(:, 2));
    % Distance can't be greater than 180 degrees, so take the min difference
    diffC(diffC>=pi) = 2*pi - diffC(diffC >= pi);

    % Compute accuracy -- threshold to 90 deg
    iCC = find(diffC>pi/2);
    AccF(i,1) = 1-(length(iCC)/length(diffC));
    diffS(i,1) = length(diffC);

    disp('Accuracy: ' + string(1-(length(iCC)/length(diffC))));
end

MeanACC = mean(AccF);