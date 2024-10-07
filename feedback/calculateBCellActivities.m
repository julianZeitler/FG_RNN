function [B1Out,B2Out] = calculateBCellActivities(E, B1, B2, G, params)
% calculateBCellActivities Calculate the activity of B-cells based on E-
%   and G-cell activities

for ori=1:params.B.numOri
    %% B1-Activity
    % Perisomatic Input (FF)
    P1 = imfilter(E.orientation(ori).data, params.B.FF.spatial_neighborhood_exc);
    % Distal Input (FB)
    % D1 = params.B.FB.scale * B1.orientation(ori).data .* (2 * ((1./(1+exp(-(imfilter(G, params.G.RF{ori+8}))))) - params.B.FB.offset)); % Distal Input (FB)
    D1 = params.B.FB.scale * (2 * ((1./(1+exp(-(imfilter(G, params.G.RF{ori+8}))))) - params.B.FB.offset)); % Distal Input (FB)

    % Normalization
    weights = gaussianFilter1DCircular(params.B.numOri, ori, 1);
    weights = max(weights) - weights;
    for i=1:params.B.numOri
        if i == 1
            OriNorm1 = B1.orientation(i).data * weights(i) + max(weights) * B2.orientation(i).data;
        else
            OriNorm1 = OriNorm1 + B1.orientation(i).data * weights(i) + max(weights) * B2.orientation(i).data;
        end
    end
    Norm1 = params.B.exp_decay + params.B.FF.inhibition * imfilter(E.orientation(ori).data, params.B.FF.spatial_neighborhood_inh) + params.B.FB.inhibition * imfilter(B2.orientation(ori).data, params.B.FB.spatial_neighborhood) + params.B.saturation * (P1.*(1 + D1)) + params.B.FF.ori_norm * OriNorm1;

    % Assign values to output and finalize calculation
    B1Out.orientation(ori).data = params.B.FF.scale*(P1.*(1 + D1))./(Norm1);
    % B1Out.orientation(ori).data = params.B.FF.scale*(P1 + D1)./(Norm1);
    B1Out.orientation(ori).data(B1Out.orientation(ori).data<0 | isnan(B1Out.orientation(ori).data)) = 0; % set invalid data to 0


    %% B2-Activity
    % Perisomatic Input (FF)
    P2 = imfilter(E.orientation(ori).data, params.B.FF.spatial_neighborhood_exc);
    % Distal Input (FB)
    % D2 = params.B.FB.scale * B2.orientation(ori).data .* (2 * ((1./(1+exp(-(imfilter(G, params.G.RF{ori}))))) - params.B.FB.offset)); % Distal Input (FB)
    D2 = params.B.FB.scale * (2 * ((1./(1+exp(-(imfilter(G, params.G.RF{ori}))))) - params.B.FB.offset)); % Distal Input (FB)

    % Normalization
    weights = gaussianFilter1DCircular(params.B.numOri, ori, 1);
    weights = max(weights) - weights;
    for i=1:params.B.numOri
        if i == 1
            OriNorm2 = B2.orientation(i).data * weights(i) + max(weights) * B1.orientation(i).data;
        else
            OriNorm2 = OriNorm2 + B2.orientation(i).data * weights(i) + max(weights) * B1.orientation(i).data;
        end
    end
    Norm2 = params.B.exp_decay + params.B.FF.inhibition * imfilter(E.orientation(ori).data, params.B.FF.spatial_neighborhood_inh) + params.B.FB.inhibition * imfilter(B1.orientation(ori).data, params.B.FB.spatial_neighborhood) + params.B.saturation * (P1.*(1 + D1)) + params.B.FF.ori_norm * OriNorm2;
    
    % Assign values to output and finalize calculation
    B2Out.orientation(ori).data = params.B.FF.scale*(P2.*(1 + D2))./(Norm2);
    % B2Out.orientation(ori).data = params.B.FF.scale*(P2 + D2)./(Norm2);
    B2Out.orientation(ori).data(B2Out.orientation(ori).data<0 | isnan(B2Out.orientation(ori).data)) = 0;
end
end

