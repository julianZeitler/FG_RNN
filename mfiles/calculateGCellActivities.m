function [G] = calculateGCellActivities(B1, B2, params)
% Calculate G-Cells with incoming B1 and B2 activities.

G = zeros(params.num_scales, size(B1, 2), size(B1, 3));
G_inh_input = zeros(params.num_scales, size(B1,2), size(B1,3));

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
                % Excitatory
                weight_exc = -cos(G_ori - B_oris(idx_B_ori));
                G(k,:,:) = squeeze(G(k,:,:)) + weight_exc.*imfilter(squeeze(B1(idx_B_ori,:,:)), params.G.RF{k, idx_G_ori});
                G(k,:,:) = squeeze(G(k,:,:)) + weight_exc.*imfilter(squeeze(B2(idx_B_ori,:,:)), params.G.RF{k, idx_G_ori+8});

                % Inhibitory
                weight_inh = cos(G_ori+pi - B_oris(idx_B_ori));
                G_inh_input(k,:,:) = squeeze(G_inh_input(k,:,:)) + weight_inh.*imfilter(squeeze(B1(idx_B_ori,:,:)), params.G.RF{k, idx_G_ori+8});
                G_inh_input(k,:,:) = squeeze(G_inh_input(k,:,:)) + weight_inh.*imfilter(squeeze(B2(idx_B_ori,:,:)), params.G.RF{k, idx_G_ori});
            elseif ori_mask2(idx_B_ori)
                % Excitatory
                weight_exc = -cos(G_ori+pi - B_oris(idx_B_ori));
                G(k,:,:) = squeeze(G(k,:,:)) + weight_exc.*imfilter(squeeze(B1(idx_B_ori,:,:)), params.G.RF{k, idx_G_ori+8});
                G(k,:,:) = squeeze(G(k,:,:)) + weight_exc.*imfilter(squeeze(B2(idx_B_ori,:,:)), params.G.RF{k, idx_G_ori});

                % Inhibitory
                weight_inh = cos(G_ori - B_oris(idx_B_ori));
                G_inh_input(k,:,:) = squeeze(G_inh_input(k,:,:)) + weight_inh.*imfilter(squeeze(B1(idx_B_ori,:,:)), params.G.RF{k, idx_G_ori});
                G_inh_input(k,:,:) = squeeze(G_inh_input(k,:,:)) + weight_inh.*imfilter(squeeze(B2(idx_B_ori,:,:)), params.G.RF{k, idx_G_ori+8});
            end
        end
    end

    G(k,:,:) = G(k,:,:)/params.num_ori;
    G(k,:,:) = params.G.scale * squeeze(G(k,:,:))./(...
        params.G.exp_decay_space + ...
        params.G.inhibition_strength_space * imfilter(squeeze(G(k,:,:)), params.G.inhibition_neighborhood{k}) + ...
        params.G.inhibitory_input_strength * squeeze(G_inh_input(k,:,:)));
end
for k=1:params.num_scales
    inhibition = zeros(size(G, 2), size(G, 3));
    for i=1:params.num_scales
        if i == k
            continue
        end
        inhibition = inhibition + imfilter(squeeze(G(i,:,:)), params.G.inhibition_neighborhood{3});
    end
    G(k,:,:) = squeeze(G(k,:,:))./...
        (params.G.exp_decay_scale + ...
        params.G.inhibition_strength_scale * inhibition);
end
end

