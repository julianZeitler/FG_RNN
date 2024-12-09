function [out] = gaussianFilter1DCircular(size, center, sigma)
% Create a gaussian-shaped filter, that has its mean at center. The
% returned array is "circular".

% Calculate values
values = zeros(size, 0);
for i=1:size
    values(i) = 1/sqrt(2*pi*sigma^2) * exp(-((i-ceil(size/2))^2)/(2*sigma^2));
end

% Assign values to corresponding index.
out = zeros(size, 0);
if center < ceil(size/2)
    for i=0:size-1 % Begin for-loop at 0 to correctly calculate mod(). Else the result for mod(size) would be 0, which is not desirable
        idx = mod(i+ceil(size/2)-center, size)+1;
        out(i+1) = values(idx);
    end
elseif center > ceil(size/2)
    for i=0:size-1
        idx = mod(i-center-ceil(size/2), size)+1;
        out(i+1) = values(idx);
    end
else
    out = values;
end
end

