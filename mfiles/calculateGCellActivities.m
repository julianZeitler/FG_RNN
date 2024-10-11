function [G] = calculateGCellActivities(B1, B2, params)
% Calculate G-Cells with incoming B1 and B2 activities.
G = zeros(params.G.num_scales, size(B1, 2), size(B2, 3));

for k=1:params.G.num_scales
    for ori=1:params.num_ori
        % collect excitatory B-cell activity from B1 and B2
        if ori==1
            G(k,:,:) = params.G.scale * (imfilter(squeeze(B1(ori,:,:)), params.G.RF{k, ori}) + imfilter(squeeze(B2(ori,:,:)), params.G.RF{k, ori+8}));
        else
            G(k,:,:) = squeeze(G(k,:,:)) + params.G.scale * (imfilter(squeeze(B1(ori,:,:)), params.G.RF{k, ori}) + imfilter(squeeze(B2(ori,:,:)), params.G.RF{k, ori+8}));
        end
        G(k,:,:) = max(0, G(k,:,:)); % remove negative values
        G(isnan(G)) = 0; % set invalid values to zero
    end
    G(k,:,:) = G(k,:,:)/params.num_ori;
    G(k,:,:) = squeeze(G(k,:,:))./(params.G.exp_decay + params.G.inhibition_strength*imfilter(squeeze(G(k,:,:)), params.G.inhibition_neighborhood)); % normalize
end
end

