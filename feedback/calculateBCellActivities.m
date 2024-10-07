function [B1,B2] = calculateBCellActivities(E, B1, B2, G, params)
% calculateBCellActivities Calculate the activity of B-cells based on E-
%   and G-cell activities

for ori=1:params.B.numOri   
    P1 = imfilter(E.orientation(ori).data, params.B.FF.spatial_neighborhood_exc); % Perisomatic Input (FF)
    % D1 = params.B.FB.scale * B1.orientation(ori).data .* (2 * ((1./(1+exp(-(imfilter(G, params.G.RF{ori+8}))))) - params.B.FB.offset)); % Distal Input (FB)
    D1 = params.B.FB.scale * (2 * ((1./(1+exp(-(imfilter(G, params.G.RF{ori+8}))))) - params.B.FB.offset)); % Distal Input (FB)
    Norm1 = params.B.exp_decay + params.B.FF.inhibition * imfilter(E.orientation(ori).data, params.B.FF.spatial_neighborhood_inh) + params.B.FB.inhibition * imfilter(B2.orientation(ori).data, params.B.FB.spatial_neighborhood) + params.B.saturation * (P1.*(1 + D1));

    P2 = imfilter(E.orientation(ori).data, params.B.FF.spatial_neighborhood_exc); % Perisomatic Input
    % D2 = params.B.FB.scale * B2.orientation(ori).data .* (2 * ((1./(1+exp(-(imfilter(G, params.G.RF{ori}))))) - params.B.FB.offset)); % Distal Input (FB)
    D2 = params.B.FB.scale * (2 * ((1./(1+exp(-(imfilter(G, params.G.RF{ori}))))) - params.B.FB.offset)); % Distal Input (FB)
    Norm2 = params.B.exp_decay + params.B.FF.inhibition * imfilter(E.orientation(ori).data, params.B.FF.spatial_neighborhood_inh) + params.B.FB.inhibition * imfilter(B1.orientation(ori).data, params.B.FB.spatial_neighborhood) + params.B.saturation * (P1.*(1 + D1));

    B1.orientation(ori).data = params.B.FF.scale*(P1.*(1 + D1))./(Norm1);
    % B1.orientation(ori).data = params.B.FF.scale*(P1 + D1)./(Norm1);
    B1.orientation(ori).data(B1.orientation(ori).data<0 | isnan(B1.orientation(ori).data)) = 0; % set invalid data to 0
    
    B2.orientation(ori).data = params.B.FF.scale*(P2.*(1 + D2))./(Norm2);
    % B2.orientation(ori).data = params.B.FF.scale*(P2 + D2)./(Norm2);
    B2.orientation(ori).data(B2.orientation(ori).data<0 | isnan(B2.orientation(ori).data)) = 0;
end
end

