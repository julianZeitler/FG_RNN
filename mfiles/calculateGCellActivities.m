function [G] = calculateGCellActivities(B1, B2, params)
% Calculate G-Cells with incoming B1 and B2 activities.
G = zeros(params.num_scales, size(B1, 2), size(B2, 3));

for k=1:params.num_scales
    for idx_G_ori = 1:params.num_ori
        G_ori = wrapTo2Pi(params.oris(idx_G_ori) + pi/2);

        % Determine the relevant orientations. Only oris +- 90° are considered
        % G_ori+pi is preferred orientation by G-cell
        ori_mask1 = abs(wrapToPi(params.oris - pi/2) - wrapToPi(G_ori+pi))<pi/2; % For oris in params.ori
        ori_mask2 = ~ori_mask1; % For opposite oris

        B_oris = wrapTo2Pi(params.oris - pi/2);

        for idx_B_ori = 1:length(B_oris)
            if ori_mask1(idx_B_ori)
                weight = -cos(G_ori - B_oris(idx_B_ori));
                G(k,:,:) = squeeze(G(k,:,:)) + weight.*imfilter(squeeze(B1(idx_B_ori,:,:)), params.G.RF{k, idx_G_ori});
                G(k,:,:) = squeeze(G(k,:,:)) + weight.*imfilter(squeeze(B2(idx_B_ori,:,:)), params.G.RF{k, idx_G_ori+8});
            elseif ori_mask2(idx_B_ori)
                weight = -cos(G_ori+pi - B_oris(idx_B_ori));
                G(k,:,:) = squeeze(G(k,:,:)) + weight.*imfilter(squeeze(B1(idx_B_ori,:,:)), params.G.RF{k, idx_G_ori+8});
                G(k,:,:) = squeeze(G(k,:,:)) + weight.*imfilter(squeeze(B2(idx_B_ori,:,:)), params.G.RF{k, idx_G_ori});
            end
        end
    end

    G(k,:,:) = G(k,:,:)/params.num_ori;
    G(k,:,:) = params.G.scale * squeeze(G(k,:,:))./...
        (params.G.exp_decay + ...
        params.G.inhibition_strength * imfilter(squeeze(G(k,:,:)), params.G.inhibition_neighborhood{k})); % normalize
end
end

