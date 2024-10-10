function [G] = calculateGCellActivities(B1, B2, params)
% Calculate G-Cells with incoming B1 and B2 activities.

for ori=1:params.numOri
    % collect excitatory B-cell activity from B1 and B2
    if ori==1
        G = params.G.scale * (imfilter(squeeze(B1(ori,:,:)), params.G.RF{ori}) + imfilter(squeeze(B2(ori,:,:)), params.G.RF{ori+8}));
    else
        G = G + params.G.scale * (imfilter(squeeze(B1(ori,:,:)), params.G.RF{ori}) + imfilter(squeeze(B2(ori,:,:)), params.G.RF{ori+8}));
    end
    % end
    G(G<0 | isnan(G)) = 0; % set invalid data to 0
end
G = G/params.numOri;
G = G./(params.G.exp_decay + params.G.inhibition_strength*imfilter(G, params.G.inhibition_neighborhood)); % normalize
end

