function [G] = calculateGCellActivities(B1, B2, params)

for ori=1:params.B.numOri
    % collect excitatory B-cell activity from B1 and B2
    if ori==1
        G = params.G.scale * (imfilter(B1.orientation(ori).data, params.G.RF{ori}) + imfilter(B2.orientation(ori).data, params.G.RF{ori+8}));
    else
        G = G + params.G.scale * (imfilter(B1.orientation(ori).data, params.G.RF{ori}) + imfilter(B2.orientation(ori).data, params.G.RF{ori+8}));
    end
    % end
    G(G<0 | isnan(G)) = 0; % set invalid data to 0
end
G = G/params.B.numOri; % normalize
% G = 10*G/max(G, [], "all");
% G = 5*2*((1./(1+exp(-G))) - 0.5);
end

