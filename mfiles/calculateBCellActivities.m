function [B1Out,B2Out] = calculateBCellActivities(E, B1, B2, G, params)
% calculateBCellActivities Calculate the activity of B-cells based on E-
%   and G-cell activities
%
% B1 and B2 are opposing B-cells. 1=Theta> and 2 = Theta<, according to
% paper notation.
B1Out = zeros(params.G.num_scales, params.num_ori, size(B1,3), size(B1,4));
B2Out = zeros(params.G.num_scales, params.num_ori, size(B2,3), size(B2,4));

for ori=1:params.num_ori
    %% B1-Activity
    % Perisomatic Input (FF)
    P1 = imfilter(squeeze(E(ori,:,:)), params.B.FF.spatial_neighborhood_exc);
    % Distal Input (FB)
    % D1 = params.B.FB.scale * B1.orientation(ori).data .* (2 * ((1./(1+exp(-(imfilter(G, params.G.RF{ori+8}))))) - params.B.FB.offset)); % Distal Input (FB)
    FB1 = zeros(params.G.num_scales, size(B1,3), size(B1,4));
    D1 = zeros(params.G.num_scales, size(B1,3), size(B1,4));
    Norm1 = zeros(params.G.num_scales, size(B1,3), size(B1,4));
    for k=1:params.G.num_scales
        FB1(k,:,:) = imfilter(squeeze(G(k,:,:)), params.G.RF{k,ori+8});
        D1(k,:,:) = params.B.FB.scale * (2 * ((1./(1+exp(-FB1(k,:,:)))) - params.B.FB.offset)); % Distal Input (FB)
    
        % FB1 = FB1/params.G.num_scales;
        
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
        Norm1(k,:,:) = params.B.exp_decay + params.B.FF.inhibition * imfilter(squeeze(E(ori,:,:)), params.B.FF.spatial_neighborhood_inh) + params.B.FB.inhibition * imfilter(squeeze(B2(k,ori,:,:)), params.B.FB.spatial_neighborhood) + params.B.saturation * (P1.*(1 + squeeze(D1(k,:,:)))) + params.B.FF.ori_norm * OriNorm1;
        
        % Assign values to output and finalize calculation
        B1Out(k,ori,:,:) = params.B.FF.scale*(P1.*(1 + squeeze(D1(k,:,:))))./(squeeze(Norm1(k,:,:)));
        % B1Out.orientation(ori).data = params.B.FF.scale*(P1 + D1)./(Norm1);
        B1Out(k,ori,:,:) = max(0, B1Out(k,ori,:,:)); % remove negative values
        B1Out(isnan(B1Out)) = 0;
    end

    %% B2-Activity
    % Perisomatic Input (FF)
    P2 = imfilter(squeeze(E(ori,:,:)), params.B.FF.spatial_neighborhood_exc);
    % Distal Input (FB)
    % D2 = params.B.FB.scale * B2.orientation(ori).data .* (2 * ((1./(1+exp(-(imfilter(G, params.G.RF{ori}))))) - params.B.FB.offset)); % Distal Input (FB)
    FB2 = zeros(params.G.num_scales, size(B2,3), size(B2,4));
    D2 = zeros(params.G.num_scales, size(B2,3), size(B2,4));
    Norm2 = zeros(params.G.num_scales, size(B2,3), size(B2,4));
    for k=1:params.G.num_scales
        FB2(k,:,:) = imfilter(squeeze(G(k,:,:)), params.G.RF{k,ori});
        D2(k,:,:) = params.B.FB.scale * (2 * ((1./(1+exp(-FB2(k,:,:)))) - params.B.FB.offset)); % Distal Input (FB)
        % FB2 = FB2/params.G.num_scales;
    
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
        Norm2(k,:,:) = params.B.exp_decay + params.B.FF.inhibition * imfilter(squeeze(E(ori,:,:)), params.B.FF.spatial_neighborhood_inh) + params.B.FB.inhibition * imfilter(squeeze(B1(k,ori,:,:)), params.B.FB.spatial_neighborhood) + params.B.saturation * (P2.*(1 + squeeze(D2(k,:,:)))) + params.B.FF.ori_norm * OriNorm2;

        % Assign values to output and finalize calculation
        B2Out(k,ori,:,:) = params.B.FF.scale*(P2.*(1 + squeeze(D2(k,:,:))))./(squeeze(Norm2(k,:,:)));
        % B2Out.orientation(ori).data = params.B.FF.scale*(P2 + D2)./(Norm2);
        B2Out(k,ori,:,:) = max(0, B2Out(k,ori,:,:)); % remove negative values
        B2Out(isnan(B2Out)) = 0;
    end
end
end

