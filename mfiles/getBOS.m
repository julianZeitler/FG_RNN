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
        % for each pixel, sum the BOS-vectors represented by their
        % orientation (edge ori +/- pi/2) and length (BOS-strength)
        % Rotate oris -90° to go from edge orientation to BOS-orientation.
        [X, Y] = pol2cart(params.oris-pi/2, (B1(:,i,j)-B2(:,i,j)).');
        [theta, rho] = cart2pol(sum(X), sum(Y)); % Perform addition in cartesian coordinates
        theta = wrapTo2Pi(theta);

        BOS(i,j,1) = theta/(2*pi); % Normalize to [0,1]
        BOS(i,j,2) = rho;
        BOS(i,j,3) = 1;

    end
end

if norm_flag==true
    BOS(:,:,2) = BOS(:,:,2)/max(BOS(:,:,2), [], "all"); % Normalize BOS-Strength
end
end

