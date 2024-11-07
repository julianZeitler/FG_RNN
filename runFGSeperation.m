function [BOS, edge_map, occ_map, group_map] = runFGSeperation(stimulus)
%RUNFGSEPERATION Summary of this function goes here
%   Detailed explanation goes here

params = makeParams();
dimensions = size(stimulus);

% Initialize activities to zero
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

for idx=1:params.iterations % Loop through all iterations
    [B1, B2] = calculateBCellActivities(E, B1, B2, G, params);
    G = calculateGCellActivities(B1, B2, params);
end

BOS = getBOS(squeeze(sum(B1,1)), squeeze(sum(B2,1)), params);
edge_map = corfresponse;
occ_map = hsv2rgb(BOS);
group_map = squeeze(sum(G));
end

