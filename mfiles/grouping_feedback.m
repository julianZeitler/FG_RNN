function [bPyr1_1, bPyr1_2, bPyr2_1, bPyr2_2, cPyr] = grouping_feedback(gPyr1,gPyr2,bPyr1_1,bPyr1_2,bPyr2_1,bPyr2_2,cPyr,params)
%Computes feedback from grouping pyramids to border-ownership cells
%
%inputs:
%gPyr - grouping cell Pyramid
%bPyr1 - Border ownership pyramid from gPyr and von Mises msk1
%bPyr2 - Border ownership pyramid from gPyr and von Mises msk2
%cPyr - complex edge Pyramid
%params - algroithm parameter structure
%
%outputs:
%bPyr1 - Border ownership pyramid from gPyr and von Mises msk1
%bPyr2 - Border ownership pyramid from gPyr and von Mises msk2
%cPyr - complex edge Pyramid
%
%By Brian Hu, Johns Hopkins University, 2017

bPrs = params.bPrs;
vmPrs = params.vmPrs;

% Similar to Alex Russell's original csPyr implementation
[vmL1, ~, vmL2, ~] = vonMisesSum(gPyr1,vmPrs); % light on dark
[vmD1, ~, vmD2, ~] = vonMisesSum(gPyr2,vmPrs); % dark on light

%calculate border ownership and grouping responses using a sigmoid function
for l = bPrs.minLevel:bPrs.maxLevel        
    for  ori = 1:bPrs.numOri

        % light on dark (similar for color channels)
        % bPyr1_1(l).orientation(ori).data = cPyr(l).orientation(ori).data .* (1+bPrs.alpha.*(1./(1+exp(-((1/1)*(vmL1(l).orientation(ori).data-vmD2(l).orientation(ori).data))))));
        P1_1 = imfilter(cPyr(l).orientation(ori).data, params.bPRs.FF.spatial_neighborhood); % Perisomatic Input (FF)
        D1_1 = params.bPRs.FB.scale * bPyr1_1(l).orientation(ori).data .* (1./(1+exp(-(vmL1(l).orientation(ori).data - params.bPRs.FB.offset)))); % Distal Input (FB)
        Norm1_1 = params.bPrs.exp_decay + params.bPRs.FF.inhibition * imfilter(cPyr(l).orientation(ori).data, params.bPRs.FB.spatial_neighborhood) + params.bPRs.FB.inhibition * imfilter(bPyr1_1(l).orientation(ori).data, params.bPRs.FB.spatial_neighborhood); % Denominator
        bPyr1_1(l).orientation(ori).data = (P1_1 + D1_1)./Norm1_1;
        bPyr1_1(l).orientation(ori).data(bPyr1_1(l).orientation(ori).data<0 | isnan(bPyr1_1(l).orientation(ori).data)) = 0; % set invalid data to 0

        % bPyr1_2(l).orientation(ori).data = cPyr(l).orientation(ori+8).data .* (1+bPrs.alpha.*(1./(1+exp(-((1/1)*(vmL2(l).orientation(ori).data-vmD1(l).orientation(ori).data))))));
        P1_2 = imfilter(cPyr(l).orientation(ori+8).data, params.bPRs.FF.spatial_neighborhood); % Perisomatic Input
        D1_2 = params.bPRs.FB.scale * bPyr1_2(l).orientation(ori).data .* (1./(1+exp(-(vmL2(l).orientation(ori).data - params.bPRs.FB.offset)))); % Distal Input
        Norm1_2 = params.bPrs.exp_decay + params.bPRs.FF.inhibition * imfilter(cPyr(l).orientation(ori+8).data, params.bPRs.FB.spatial_neighborhood) + params.bPRs.FB.inhibition * imfilter(bPyr1_2(l).orientation(ori).data, params.bPRs.FB.spatial_neighborhood); % Denominator
        bPyr1_2(l).orientation(ori).data = (P1_2 + D1_2)./Norm1_2;
        bPyr1_2(l).orientation(ori).data(bPyr1_2(l).orientation(ori).data<0 | isnan(bPyr1_2(l).orientation(ori).data)) = 0;

        % dark on light (similar for color channels)
        % bPyr2_1(l).orientation(ori).data = cPyr(l).orientation(ori+8).data .* (1+bPrs.alpha.*(1./(1+exp(-((1/1)*(vmD1(l).orientation(ori).data-vmL2(l).orientation(ori).data))))));
        P2_1 = imfilter(cPyr(l).orientation(ori+8).data, params.bPRs.FF.spatial_neighborhood); % Perisomatic Input
        D2_1 = params.bPRs.FB.scale * bPyr2_1(l).orientation(ori).data .* (1./(1+exp(-(vmD1(l).orientation(ori).data - params.bPRs.FB.offset)))); % Distal Input
        Norm2_1 = params.bPrs.exp_decay + params.bPRs.FF.inhibition * imfilter(cPyr(l).orientation(ori+8).data, params.bPRs.FB.spatial_neighborhood) + params.bPRs.FB.inhibition * imfilter(bPyr2_1(l).orientation(ori).data, params.bPRs.FB.spatial_neighborhood); % Denominator
        bPyr2_1(l).orientation(ori).data = (P2_1 + D2_1)./Norm2_1;
        bPyr2_1(l).orientation(ori).data(bPyr2_1(l).orientation(ori).data<0 | isnan(bPyr2_1(l).orientation(ori).data)) = 0;

        % bPyr2_2(l).orientation(ori).data = cPyr(l).orientation(ori).data .* (1+bPrs.alpha.*(1./(1+exp(-((1/1)*(vmD2(l).orientation(ori).data-vmL1(l).orientation(ori).data))))));
        P2_2 = imfilter(cPyr(l).orientation(ori).data, params.bPRs.FF.spatial_neighborhood); % Perisomatic Input
        D2_2 = params.bPRs.FB.scale * bPyr2_2(l).orientation(ori).data .* (1./(1+exp(-(vmD2(l).orientation(ori).data - params.bPRs.FB.offset)))); % Distal Input
        Norm2_2 = params.bPrs.exp_decay + params.bPRs.FF.inhibition * imfilter(cPyr(l).orientation(ori).data, params.bPRs.FB.spatial_neighborhood) + params.bPRs.FB.inhibition * imfilter(bPyr2_2(l).orientation(ori).data, params.bPRs.FB.spatial_neighborhood); % Denominator
        bPyr2_2(l).orientation(ori).data = (P2_2 + D2_2)./Norm2_2;
        bPyr2_2(l).orientation(ori).data(bPyr2_2(l).orientation(ori).data<0 | isnan(bPyr2_2(l).orientation(ori).data)) = 0;

    end
end

