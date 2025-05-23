function [B1Out,B2Out] = calculateBCellActivities(E, B1, B2, G, params, debug)
% calculateBCellActivities Calculate the activity of B-cells based on E-
%   and G-cell activities
%
% B1 and B2 are opposing B-cells. 1=Theta> and 2 = Theta<, according to
% paper notation.
B1Out = zeros(params.num_ori, size(B1,2), size(B1,3));
B2Out = zeros(params.num_ori, size(B2,2), size(B2,3));

for ori=1:params.num_ori
    % Perisomatic Input (FF)
    P = imfilter(squeeze(E(ori,:,:)), params.B.FF.spatial_neighborhood_exc);

    %% B1-Activity
    % Calculate FB-Signal
    FB1 = zeros(params.num_scales, size(B1,2), size(B1,3));
    for k=1:params.num_scales
        B_ori = wrapToPi(params.oris(ori)-pi/2);
        G_oris = wrapToPi(params.oris + pi/2);

        % Determine the relevant orientations. Only oris +- 90° are considered
        % B_ori is preferred orientation by G-cell
        ori_mask1 = abs(wrapToPi(G_oris - B_ori+pi))<pi/2; % For oris in params.ori
        ori_mask2 = ~ori_mask1; % For opposite oris        

        for idx_G_ori = 1:length(G_oris)
            if ori_mask1(idx_G_ori)
                weight = -cos(B_ori - G_oris(idx_G_ori));
                FB1(k,:,:) = squeeze(FB1(k,:,:)) + weight.*imfilter(squeeze(G(k,:,:)), params.G.RF{k, idx_G_ori+8});
            elseif ori_mask2(idx_G_ori)
                weight = -cos(B_ori+pi - G_oris(idx_G_ori));
                FB1(k,:,:) = squeeze(FB1(k,:,:)) + weight.*imfilter(squeeze(G(k,:,:)), params.G.RF{k, idx_G_ori});
            end
        end
    end
    
    [FB1, FB1_idx] = max(FB1, [], 1);
    FB1 = squeeze(FB1);

    % Distal Input (FB)
    D1 = params.B.FB.scale * (2 * ((1./(1+exp(-FB1))) - params.B.FB.offset)); % Distal Input (FB)

    % Normalization
    weights = gaussianFilter1DCircular(params.num_ori, ori, 1);
    weights = max(weights) - weights;
    borderSupression = zeros(size(B1,2), size(B1,3));
    for i=1:params.num_ori
        if i == 1
            OriNorm1 = squeeze(B1(i,:,:)) * weights(i) + max(weights) * squeeze(B2(i,:,:));
        else
            OriNorm1 = OriNorm1 + squeeze(B1(i,:,:)) * weights(i) + max(weights) * squeeze(B2(i,:,:));
        end
        borderSupression = borderSupression + imfilter(squeeze(B1(i,:,:)), borderSupressionKernel(20, params.oris(i)));
        borderSupression = borderSupression + imfilter(squeeze(B2(i,:,:)), borderSupressionKernel(20, params.oris(i)+pi));
    end
    
    Norm1 = params.B.exp_decay + ...
        params.B.FF.inhibition * imfilter(squeeze(E(ori,:,:)), params.B.FF.spatial_neighborhood_inh) + ...
        params.B.FB.inhibition * imfilter(squeeze(B2(ori,:,:)), params.B.FB.spatial_neighborhood) + ...
        params.B.saturation * (P.*(1 + D1)) + ...
        params.B.FF.ori_norm * OriNorm1 + ...
        params.B.FB.border_supression * borderSupression;
    
    % Assign values to output and finalize calculation
    B1Out(ori,:,:) = params.B.FF.scale*(P.*(1 + D1))./Norm1;
    B1Out(ori,:,:) = max(0, B1Out(ori,:,:)); % remove negative values
    B1Out(isnan(B1Out)) = 0;

    %% B2-Activity
    % Calculate FB-Signal
    FB2 = zeros(params.num_scales, size(B2,2), size(B2,3));
    for k=1:params.num_scales
        B_ori = wrapToPi(params.oris(ori)-pi/2);
        G_oris = wrapToPi(params.oris + pi/2);

        % Determine the relevant orientations. Only oris +- 90° are considered
        % B_ori is preferred orientation by G-cell
        ori_mask1 = abs(wrapToPi(G_oris - B_ori+pi))<pi/2; % For oris in params.ori
        ori_mask2 = ~ori_mask1; % For opposite oris

        for idx_G_ori = 1:length(G_oris)
            if ori_mask1(idx_G_ori)
                weight = -cos(B_ori - G_oris(idx_G_ori));
                FB2(k,:,:) = squeeze(FB2(k,:,:)) + weight.*imfilter(squeeze(G(k,:,:)), params.G.RF{k, idx_G_ori});
            elseif ori_mask2(idx_G_ori)
                weight = -cos(B_ori+pi - G_oris(idx_G_ori));
                FB2(k,:,:) = squeeze(FB2(k,:,:)) + weight.*imfilter(squeeze(G(k,:,:)), params.G.RF{k, idx_G_ori+8});
            end
        end
    end
    
    [FB2, FB2_idx] = max(FB2, [], 1);
    FB2 = squeeze(FB2);

    % Distal Input (FB)
    D2 = params.B.FB.scale * (2 * ((1./(1+exp(-FB2))) - params.B.FB.offset)); % Distal Input (FB)

    % Normalization
    weights = gaussianFilter1DCircular(params.num_ori, ori, 1);
    weights = max(weights) - weights;
    borderSupression = zeros(size(B1,2), size(B1,3));
    for i=1:params.num_ori
        if i == 1
            OriNorm2 = squeeze(B2(i,:,:)) * weights(i) + max(weights) * squeeze(B1(i,:,:));
        else
            OriNorm2 = OriNorm2 + squeeze(B2(i,:,:)) * weights(i) + max(weights) * squeeze(B1(i,:,:));
        end
        borderSupression = borderSupression + imfilter(squeeze(B1(i,:,:)), borderSupressionKernel(20, params.oris(i)+pi));
        borderSupression = borderSupression + imfilter(squeeze(B2(i,:,:)), borderSupressionKernel(20, params.oris(i)));
    end

    Norm2 = params.B.exp_decay + ...
        params.B.FF.inhibition * imfilter(squeeze(E(ori,:,:)), params.B.FF.spatial_neighborhood_inh) + ...
        params.B.FB.inhibition * imfilter(squeeze(B1(ori,:,:)), params.B.FB.spatial_neighborhood) + ...
        params.B.saturation * (P.*(1 + D2)) + ...
        params.B.FF.ori_norm * OriNorm2 + ...
        params.B.FB.border_supression * borderSupression;

    % Assign values to output and finalize calculation
    B2Out(ori,:,:) = params.B.FF.scale*(P.*(1 + D2))./Norm2;
    B2Out(ori,:,:) = max(0, B2Out(ori,:,:)); % remove negative values
    B2Out(isnan(B2Out)) = 0;

    %% Debug
    if debug == true
        % Compare B1 and B2
        B1_vs_B2_mask = squeeze(B1Out(ori,:,:) > B2Out(ori,:,:)); % 1 if B1 wins, 0 if B2 wins
        
        % Store the best B1 or B2 activity
        BestActivityAtOri = squeeze(B1Out(ori,:,:));
        B2Map = squeeze(B2Out(ori,:,:)); % 2D matrix
        BestActivityAtOri(~B1_vs_B2_mask) = B2Map(~B1_vs_B2_mask);
        
        % Store scale indices
        FB1_idx = squeeze(FB1_idx);
        FB2_idx = squeeze(FB2_idx);
        BestScaleIdxAtOri = FB1_idx; % default
        BestScaleIdxAtOri(~B1_vs_B2_mask) = FB2_idx(~B1_vs_B2_mask);
        
        % Save per orientation
        BestActivity(:,:,:,ori) = BestActivityAtOri;
        BestScaleIdx(:,:,:,ori) = BestScaleIdxAtOri;
    end
end

if debug == true
    % Find maximum over orientations
    [~, WinningOri] = max(BestActivity, [], 4); % 4th dimension = orientation
    
    % Initialize final map
    [H, W, ~, ~] = size(BestActivity);
    FinalWinningScaleIdx = zeros(H, W);
    
    % Go through each pixel and select the winning scale
    for x = 1:H
        for y = 1:W
            FinalWinningScaleIdx(x,y) = BestScaleIdx(x,y,1,WinningOri(x,y));
        end
    end
    
    % Plot
    figure;
    imagesc(FinalWinningScaleIdx);
    axis image;
    colorbar;
    colormap(jet(params.num_scales)); % or your preferred map
    title('Final Winning Scale Index');
    xlabel('X');
    ylabel('Y');
end
end

