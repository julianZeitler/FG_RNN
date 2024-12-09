function [contrastStimulus] = verticalContrastStimulus(dim1, dim2)
%VERTICALCONTRAST Summary of this function goes here
%   Detailed explanation goes here
contrastStimulus = zeros(dim1, dim2); 
contrastStimulus(:, (dim2/2+1):end) = 1;
end

