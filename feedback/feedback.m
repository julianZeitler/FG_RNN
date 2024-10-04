function [BL1, BL2, BD1, BD2] = feedback(GL, GD, BL1, BL2, BD1, BD2, E, params)
%Computes feedback from grouping cells to border-ownership cells

%calculate B-cell activity based on two-point integration scheme
for  ori = 1:params.B.numOri
    % light on dark (similar for color channels)
    P1_1 = imfilter(E.orientation(ori).data, params.B.FF.spatial_neighborhood); % Perisomatic Input (FF)
    D1_1 = params.B.FB.scale * BL1.orientation(ori).data .* (1./(1+exp(-(imfilter(GL, params.vmPrs.msk_1{ori}) - params.B.FB.offset)))); % Distal Input (FB)
    Norm1_1 = params.B.exp_decay + params.B.FF.inhibition * imfilter(E.orientation(ori).data, params.B.FB.spatial_neighborhood) + params.B.FB.inhibition * imfilter(BL2.orientation(ori).data, params.B.FB.spatial_neighborhood); % Denominator
    BL1.orientation(ori).data = (P1_1 + D1_1)./Norm1_1;
    BL1.orientation(ori).data(BL1.orientation(ori).data<0 | isnan(BL1.orientation(ori).data)) = 0; % set invalid data to 0

    P1_2 = imfilter(E.orientation(ori+8).data, params.B.FF.spatial_neighborhood); % Perisomatic Input
    D1_2 = params.B.FB.scale * BL2.orientation(ori).data .* (1./(1+exp(-(imfilter(GL, params.vmPrs.msk_2{ori}) - params.B.FB.offset)))); % Distal Input
    Norm1_2 = params.B.exp_decay + params.B.FF.inhibition * imfilter(E.orientation(ori+8).data, params.B.FB.spatial_neighborhood) + params.B.FB.inhibition * imfilter(BL1.orientation(ori).data, params.B.FB.spatial_neighborhood); % Denominator
    BL2.orientation(ori).data = (P1_2 + D1_2)./Norm1_2;
    BL2.orientation(ori).data(BL2.orientation(ori).data<0 | isnan(BL2.orientation(ori).data)) = 0;

    % dark on light (similar for color channels)
    P2_1 = imfilter(E.orientation(ori+8).data, params.B.FF.spatial_neighborhood); % Perisomatic Input
    D2_1 = params.B.FB.scale * BD1(l).orientation(ori).data .* (1./(1+exp(-(imfilter(GD, params.vmPrs.msk_1{ori}) - params.B.FB.offset)))); % Distal Input
    Norm2_1 = params.B.exp_decay + params.B.FF.inhibition * imfilter(E.orientation(ori+8).data, params.B.FB.spatial_neighborhood) + params.B.FB.inhibition * imfilter(BD2(l).orientation(ori).data, params.B.FB.spatial_neighborhood); % Denominator
    BD1(l).orientation(ori).data = (P2_1 + D2_1)./Norm2_1;
    BD1(l).orientation(ori).data(BD1(l).orientation(ori).data<0 | isnan(BD1(l).orientation(ori).data)) = 0;

    P2_2 = imfilter(E.orientation(ori).data, params.B.FF.spatial_neighborhood); % Perisomatic Input
    D2_2 = params.B.FB.scale * BD2(l).orientation(ori).data .* (1./(1+exp(-(imfilter(GD, params.vmPrs.msk_2{ori}) - params.B.FB.offset)))); % Distal Input
    Norm2_2 = params.B.exp_decay + params.B.FF.inhibition * imfilter(E.orientation(ori).data, params.B.FB.spatial_neighborhood) + params.B.FB.inhibition * imfilter(BD1(l).orientation(ori).data, params.B.FB.spatial_neighborhood); % Denominator
    BD2(l).orientation(ori).data = (P2_2 + D2_2)./Norm2_2;
    BD2(l).orientation(ori).data(BD2(l).orientation(ori).data<0 | isnan(BD2(l).orientation(ori).data)) = 0;

    % plot feedback activity
    % if l == 1
    %     figure;
    %     t = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    %     sgtitle(['Orientation ', num2str(ori)]);
    % 
    %     nexttile; 
    %     imagesc(D1_1); % Plot D1_1
    %     title(['FB L Ori ', num2str(ori)]);
    %     axis off; colorbar;
    % 
    %     nexttile; 
    %     imagesc(D1_2); % Plot D1_2
    %     title(['FB L Ori ', num2str(ori), ' + pi']);
    %     axis off; colorbar;
    % 
    %     nexttile; 
    %     imagesc(D2_1); % Plot D2_1
    %     title(['FB D Ori ', num2str(ori)]);
    %     axis off; colorbar;
    % 
    %     nexttile; 
    %     imagesc(D2_2); % Plot D2_2
    %     title(['FB D Ori ', num2str(ori), ' + pi']);
    %     axis off; colorbar
    % end
end

