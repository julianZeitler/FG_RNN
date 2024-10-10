function [BOS] = getBOS(B1, B2, params, norm)
%COMBINEBOS Calculate and combine the BOS-Signal
%   Inputs:
%       B1: oxmxn Matrix (o=#oris, mxn=model dimensions)
%       B2: oxmxn Matrix (o=#oris, mxn=model dimensions)
%       params: struct created by "makeParams.m"
%       norm: optional ["normalized", "unnormalized"], default is
%       "normalized"
%
%   Out:
%       BOS: mxnx3 Matrix (mxn=model dimensions, 3: HSV-values)
%           H - range [0,1] (encodes [0,2pi]) -> BOS-direction
%           S - normalized strength of BOS-Signal
%           V - 1
norm_flag = true;
if nargin > 3
    if norm=="normalized"
        norm_flag = true;
    elseif norm=="unnormalized"
        norm_flag = false;
    end
end

dimensions = [size(B1, 2), size(B2, 3)];
BOS = zeros(dimensions(1), dimensions(2), 3);

for i=1:dimensions(1)
    for j=1:dimensions(2)
        % for each pixel, get the maximum BOS-Value. This direction will be
        % used to calculate the H-Value
        [max_val,max_idx] = max(B1(:,i,j)-B2(:,i,j));
        [min_val,min_idx] = min(B1(:,i,j)-B2(:,i,j));
        if abs(max_val) > abs(min_val)
            % BOS-direction is -90° to edge orientation
            theta1 = params.oris(max_idx)-pi/2;
            if theta1 < 0
                theta1 = theta1 + 2*pi;
            elseif theta1 > 2*pi
                theta1 = theta1 - 2*pi;
            end
            BOS(i,j,1) = (theta1)/(2*pi); % Normalize to [0,1]
            BOS(i,j,2) = abs(max_val);
            BOS(i,j,3) = 1;
        else
            % BOS-direction is +90° to edge orientation
            theta2 = params.oris(min_idx)+pi/2;
            if theta2 < 0
                theta2 = theta2 + 2*pi;
            elseif theta2 > 2*pi
                theta2 = theta2 - 2*pi;
            end
            BOS(i,j,1) = (theta2)/(2*pi); % Normalize to [0,1]
            BOS(i,j,2) = abs(min_val);
            BOS(i,j,3) = 1;
        end

    end
end

if norm_flag==true
    BOS(:,:,2) = BOS(:,:,2)/max(BOS(:,:,2), [], "all"); % Normalize BOS-Strength
end
end

