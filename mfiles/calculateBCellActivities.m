function [B1Out,B2Out] = calculateBCellActivities(E, B1, B2, G, params)
% calculateBCellActivities Calculate the activity of B-cells based on E-
%   and G-cell activities
%
% B1 and B2 are opposing B-cells. 1=Theta> and 2 = Theta<, according to
% paper notation.

for ori=1:params.numOri
    %% B1-Activity
    % Perisomatic Input (FF)
    P1 = imfilter(squeeze(E(ori,:,:)), params.B.FF.spatial_neighborhood_exc);
    % Distal Input (FB)
    % D1 = params.B.FB.scale * B1.orientation(ori).data .* (2 * ((1./(1+exp(-(imfilter(G, params.G.RF{ori+8}))))) - params.B.FB.offset)); % Distal Input (FB)
    D1 = params.B.FB.scale * (2 * ((1./(1+exp(-(imfilter(G, params.G.RF{ori+8}))))) - params.B.FB.offset)); % Distal Input (FB)

    % Normalization
    weights = gaussianFilter1DCircular(params.numOri, ori, 1);
    weights = max(weights) - weights;
    for i=1:params.numOri
        if i == 1
            OriNorm1 = squeeze(B1(i,:,:)) * weights(i) + max(weights) * squeeze(B2(i,:,:));
        else
            OriNorm1 = OriNorm1 + squeeze(B1(i,:,:)) * weights(i) + max(weights) * squeeze(B2(i,:,:));
        end
    end
    Norm1 = params.B.exp_decay + params.B.FF.inhibition * imfilter(squeeze(E(ori,:,:)), params.B.FF.spatial_neighborhood_inh) + params.B.FB.inhibition * imfilter(squeeze(B2(ori,:,:)), params.B.FB.spatial_neighborhood) + params.B.saturation * (P1.*(1 + D1)) + params.B.FF.ori_norm * OriNorm1;

    % Assign values to output and finalize calculation
    B1Out(ori,:,:) = params.B.FF.scale*(P1.*(1 + D1))./(Norm1);
    % B1Out.orientation(ori).data = params.B.FF.scale*(P1 + D1)./(Norm1);
    B1(ori,:,:) = max(0, B1(ori,:,:)); % remove negative values

    %% B2-Activity
    % Perisomatic Input (FF)
    P2 = imfilter(squeeze(E(ori,:,:)), params.B.FF.spatial_neighborhood_exc);
    % Distal Input (FB)
    % D2 = params.B.FB.scale * B2.orientation(ori).data .* (2 * ((1./(1+exp(-(imfilter(G, params.G.RF{ori}))))) - params.B.FB.offset)); % Distal Input (FB)
    D2 = params.B.FB.scale * (2 * ((1./(1+exp(-(imfilter(G, params.G.RF{ori}))))) - params.B.FB.offset)); % Distal Input (FB)

    % Normalization
    weights = gaussianFilter1DCircular(params.numOri, ori, 1);
    weights = max(weights) - weights;
    for i=1:params.numOri
        if i == 1
            OriNorm2 = squeeze(B2(i,:,:)) * weights(i) + max(weights) * squeeze(B1(i,:,:));
        else
            OriNorm2 = OriNorm2 + squeeze(B2(i,:,:)) * weights(i) + max(weights) * squeeze(B1(i,:,:));
        end
    end
    Norm2 = params.B.exp_decay + params.B.FF.inhibition * imfilter(squeeze(E(ori,:,:)), params.B.FF.spatial_neighborhood_inh) + params.B.FB.inhibition * imfilter(squeeze(B1(ori,:,:)), params.B.FB.spatial_neighborhood) + params.B.saturation * (P2.*(1 + D2)) + params.B.FF.ori_norm * OriNorm2;
    
    % Assign values to output and finalize calculation
    B2Out(ori,:,:) = params.B.FF.scale*(P2.*(1 + D2))./(Norm2);
    % B2Out.orientation(ori).data = params.B.FF.scale*(P2 + D2)./(Norm2);
    B2Out(ori,:,:) = max(0, B2Out(ori,:,:)); % remove negative values
end
end

