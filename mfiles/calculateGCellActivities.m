function [G] = calculateGCellActivities(B1, B2, params, attention)
% Calculate G-Cells with incoming B1 and B2 activities.

if nargin < 4
    attention = zeros(1, params.num_scales);
end

G = zeros(params.num_scales, size(B1, 2), size(B1, 3));
G_exc_input = zeros(params.num_scales, size(B1,2), size(B1,3));
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
                G_exc_input(k,:,:) = squeeze(G_exc_input(k,:,:)) + weight_exc.*imfilter(squeeze(B1(idx_B_ori,:,:)), params.G.RF{k, idx_G_ori});
                G_exc_input(k,:,:) = squeeze(G_exc_input(k,:,:)) + weight_exc.*imfilter(squeeze(B2(idx_B_ori,:,:)), params.G.RF{k, idx_G_ori+8});

                % Inhibitory
                weight_inh = cos(G_ori+pi - B_oris(idx_B_ori));
                G_inh_input(k,:,:) = squeeze(G_inh_input(k,:,:)) + weight_inh.*imfilter(squeeze(B1(idx_B_ori,:,:)), params.G.RF{k, idx_G_ori+8});
                G_inh_input(k,:,:) = squeeze(G_inh_input(k,:,:)) + weight_inh.*imfilter(squeeze(B2(idx_B_ori,:,:)), params.G.RF{k, idx_G_ori});
            elseif ori_mask2(idx_B_ori)
                % Excitatory
                weight_exc = -cos(G_ori+pi - B_oris(idx_B_ori));
                G_exc_input(k,:,:) = squeeze(G_exc_input(k,:,:)) + weight_exc.*imfilter(squeeze(B1(idx_B_ori,:,:)), params.G.RF{k, idx_G_ori+8});
                G_exc_input(k,:,:) = squeeze(G_exc_input(k,:,:)) + weight_exc.*imfilter(squeeze(B2(idx_B_ori,:,:)), params.G.RF{k, idx_G_ori});

                % Inhibitory
                weight_inh = cos(G_ori - B_oris(idx_B_ori));
                G_inh_input(k,:,:) = squeeze(G_inh_input(k,:,:)) + weight_inh.*imfilter(squeeze(B1(idx_B_ori,:,:)), params.G.RF{k, idx_G_ori});
                G_inh_input(k,:,:) = squeeze(G_inh_input(k,:,:)) + weight_inh.*imfilter(squeeze(B2(idx_B_ori,:,:)), params.G.RF{k, idx_G_ori+8});
            end
        end
    end
    G_exc_input(k,:,:) = G_exc_input(k,:,:)/params.num_ori;
    G_inh_input(k,:,:) = G_inh_input(k,:,:)/params.num_ori;
end

for k=1:params.num_scales
    scale_inhibition = zeros(size(G, 2), size(G, 3));
    for l=1:params.num_scales
        if l == k
            continue
        end
        scale_inhibition = scale_inhibition + squeeze(G_exc_input(l,:,:)); %imfilter(squeeze(G_exc_input(l,:,:)), params.G.inhibition_neighborhood{1});
    end
    G(k,:,:) = params.G.beta * squeeze(G_exc_input(k,:,:))*(1+attention(k))./(...
        params.G.alpha + ...
        params.G.gamma * imfilter(squeeze(G_exc_input(k,:,:)), params.G.inhibition_neighborhood{k}) + ...
        params.G.zeta * squeeze(G_exc_input(k,:,:)) + ...
        params.G.mu * squeeze(G_inh_input(k,:,:)) + ...
        params.G.nu * scale_inhibition);
end
end

