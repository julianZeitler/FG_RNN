function [B1Out,B2Out] = calculateBCellActivities(E, B1, B2, G, params)
% calculateBCellActivities Calculate the activity of B-cells based on E-
%   and G-cell activities
%
% B1 and B2 are opposing B-cells. 1=Theta> and 2 = Theta<, according to
% paper notation.
B1Out = zeros(params.num_scales, params.num_ori, size(B1,3), size(B1,4));
B2Out = zeros(params.num_scales, params.num_ori, size(B2,3), size(B2,4));

for ori=1:params.num_ori
    % Perisomatic Input (FF)
    P = zeros(params.num_scales, size(B1,3), size(B1,4));
    P(1,:,:) = imfilter(squeeze(E(ori,:,:)), params.B.FF.spatial_neighborhood_exc);
    for k=2:params.num_scales
        P(k,:,:) = imgaussfilt(squeeze(P(1,:,:)), params.B.FB.coarse_scale*params.scale_step^(k-1));
    end

    %% B1-Activity
    FB1 = zeros(params.num_scales, size(B1,3), size(B1,4));
    for k=1:params.num_scales
        FB1(k,:,:) = imfilter(squeeze(G(k,:,:)), params.G.RF{k,ori+8});
    end
    for k=1:params.num_scales
        for j=k+1:params.num_scales
            FB1(k,:,:) = FB1(k,:,:) + 1/(2^(j-k+1))*FB1(j,:,:);
        end
    end

    D1 = zeros(params.num_scales, size(B1,3), size(B1,4));
    Norm1 = zeros(params.num_scales, size(B1,3), size(B1,4));
    for k=1:params.num_scales
        % Distal Input (FB)
        D1(k,:,:) = params.B.FB.scale * (2 * ((1./(1+exp(-FB1(k,:,:)))) - params.B.FB.offset)); % Distal Input (FB)
    
        % FB1 = FB1/params.num_scales;
        
        % Normalization
        weights = gaussianFilter1DCircular(params.num_ori, ori, 1);
        weights = max(weights) - weights;
        for i=1:params.num_ori
            if i == 1
                OriNorm1 = squeeze(B1(k,i,:,:)) * weights(i) + max(weights) * squeeze(B2(k,i,:,:));
            else
                OriNorm1 = OriNorm1 + squeeze(B1(k,i,:,:)) * weights(i) + max(weights) * squeeze(B2(k,i,:,:));
            end
        end
        Norm1(k,:,:) = params.B.exp_decay + params.B.FF.inhibition * imfilter(squeeze(E(ori,:,:)), params.B.FF.spatial_neighborhood_inh) + params.B.FB.inhibition * imfilter(squeeze(B2(k,ori,:,:)), params.B.FB.spatial_neighborhood) + params.B.saturation * (squeeze(P(k,:,:)).*(1 + squeeze(D1(k,:,:)))) + params.B.FF.ori_norm * OriNorm1;
        
        % Assign values to output and finalize calculation
        B1Out(k,ori,:,:) = params.B.FF.scale*(squeeze(P(k,:,:)).*(1 + squeeze(D1(k,:,:))))./(squeeze(Norm1(k,:,:)));
        B1Out(k,ori,:,:) = max(0, B1Out(k,ori,:,:)); % remove negative values
        B1Out(isnan(B1Out)) = 0;
    end

    %% B2-Activity
    FB2 = zeros(params.num_scales, size(B2,3), size(B2,4));
    for k=1:params.num_scales
        FB2(k,:,:) = imfilter(squeeze(G(k,:,:)), params.G.RF{k,ori});
    end
    for k=1:params.num_scales
        for j=k+1:params.num_scales
            FB2(k,:,:) = FB2(k,:,:) + 1/(2^(j-k+1))*FB2(j,:,:);
        end
    end
    D2 = zeros(params.num_scales, size(B2,3), size(B2,4));
    Norm2 = zeros(params.num_scales, size(B2,3), size(B2,4));
    for k=1:params.num_scales
        % Distal Input (FB)
        D2(k,:,:) = params.B.FB.scale * (2 * ((1./(1+exp(-FB2(k,:,:)))) - params.B.FB.offset)); % Distal Input (FB)
        % FB2 = FB2/params.num_scales;
    
        % Normalization
        weights = gaussianFilter1DCircular(params.num_ori, ori, 1);
        weights = max(weights) - weights;
        for i=1:params.num_ori
            if i == 1
                OriNorm2 = squeeze(B2(k,i,:,:)) * weights(i) + max(weights) * squeeze(B1(k,i,:,:));
            else
                OriNorm2 = OriNorm2 + squeeze(B2(k,i,:,:)) * weights(i) + max(weights) * squeeze(B1(k,i,:,:));
            end
        end
        Norm2(k,:,:) = params.B.exp_decay + params.B.FF.inhibition * imfilter(squeeze(E(ori,:,:)), params.B.FF.spatial_neighborhood_inh) + params.B.FB.inhibition * imfilter(squeeze(B1(k,ori,:,:)), params.B.FB.spatial_neighborhood) + params.B.saturation * (squeeze(P(k,:,:)).*(1 + squeeze(D2(k,:,:)))) + params.B.FF.ori_norm * OriNorm2;

        % Assign values to output and finalize calculation
        B2Out(k,ori,:,:) = params.B.FF.scale*(squeeze(P(k,:,:)).*(1 + squeeze(D2(k,:,:))))./(squeeze(Norm2(k,:,:)));
        B2Out(k,ori,:,:) = max(0, B2Out(k,ori,:,:)); % remove negative values
        B2Out(isnan(B2Out)) = 0;
    end
end
end

